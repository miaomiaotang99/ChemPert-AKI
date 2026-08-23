"""Trainable Morgan + directed-MPNN chemical encoder (no PyG dependency)."""
from dataclasses import dataclass
from typing import Literal
import torch
from torch import nn
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

ATOM_DIM, BOND_DIM = 38, 10

def _oh(v, xs): return [float(v == x) for x in xs] + [float(v not in xs)]
def _atom(a):
    x = _oh(a.GetAtomicNum(), [1,6,7,8,9,15,16,17,35,53])
    x += _oh(a.GetDegree(), [0,1,2,3,4,5]) + _oh(a.GetFormalCharge(), [-2,-1,0,1,2])
    x += _oh(int(a.GetHybridization()), [2,3,4,5,6]) + _oh(a.GetTotalNumHs(), [0,1,2,3,4])
    return x + [float(a.GetIsAromatic()), float(a.IsInRing())]
def _bond(b):
    x = [float(b.GetBondType() == t) for t in (Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE, Chem.BondType.AROMATIC)]
    return x + [float(b.GetIsConjugated()), float(b.IsInRing())] + _oh(int(b.GetStereo()), [0,2,3])

@dataclass
class MoleculeBatch:
    morgan: torch.Tensor; atom_x: torch.Tensor; bond_x: torch.Tensor
    edge_index: torch.Tensor; reverse_edge: torch.Tensor; atom_batch: torch.Tensor
    def to(self, device):
        return MoleculeBatch(**{k: v.to(device) for k, v in vars(self).items()})

def collate_molecules(smiles, radius=2, n_bits=2048):
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    fps=[]; atoms=[]; bonds=[]; edges=[]; rev=[]; groups=[]; ao=eo=0
    for bi, smi in enumerate(smiles):
        mol=Chem.MolFromSmiles(str(smi))
        if mol is None or not mol.GetNumAtoms(): raise ValueError(f"Invalid/empty SMILES at index {bi}: {smi!r}")
        fp=torch.zeros(n_bits); fp[list(gen.GetFingerprint(mol).GetOnBits())]=1; fps.append(fp)
        atoms += [_atom(a) for a in mol.GetAtoms()]; groups += [bi]*mol.GetNumAtoms()
        for b in mol.GetBonds():
            u,v=b.GetBeginAtomIdx()+ao,b.GetEndAtomIdx()+ao; bf=_bond(b)
            edges += [(u,v),(v,u)]; bonds += [bf,bf]; rev += [eo+1,eo]; eo+=2
        ao += mol.GetNumAtoms()
    ei=torch.tensor(edges,dtype=torch.long).t().contiguous() if edges else torch.empty(2,0,dtype=torch.long)
    return MoleculeBatch(torch.stack(fps),torch.tensor(atoms),torch.tensor(bonds,dtype=torch.float32).reshape(-1,BOND_DIM),ei,torch.tensor(rev),torch.tensor(groups))

@dataclass
class ChemicalEncoderConfig:
    mode: Literal["morgan","graph","concat","cross_attention"]="cross_attention"
    n_bits:int=2048; hidden_dim:int=256; output_dim:int=256; message_steps:int=3; num_heads:int=8; dropout:float=.1

def _pad(x,g,b):
    counts=torch.bincount(g,minlength=b); n=int(counts.max()); y=x.new_zeros(b,n,x.shape[-1]); m=torch.ones(b,n,dtype=torch.bool,device=x.device); pos=torch.zeros(b,dtype=torch.long,device=x.device)
    for i in range(len(x)): j=int(g[i]); k=int(pos[j]); y[j,k]=x[i]; m[j,k]=False; pos[j]+=1
    return y,m
def _mean(x,m):
    z=(~m).unsqueeze(-1); return (x*z).sum(1)/z.sum(1).clamp_min(1)

class _Morgan(nn.Module):
    def __init__(self,c): super().__init__(); self.n=c.n_bits; self.e=nn.Embedding(c.n_bits+1,c.hidden_dim,padding_idx=c.n_bits); self.nm=nn.LayerNorm(c.hidden_dim); self.d=nn.Dropout(c.dropout)
    def forward(self,fp):
        ids=[torch.where(x>0)[0] for x in fp]; n=max(1,max(map(len,ids))); out=torch.full((len(ids),n),self.n,device=fp.device); mask=torch.ones_like(out,dtype=torch.bool)
        for i,x in enumerate(ids): out[i,:len(x)]=x; mask[i,:len(x)]=False
        return self.d(self.nm(self.e(out.long()))),mask
class _DMPNN(nn.Module):
    def __init__(self,c): super().__init__(); h=c.hidden_dim; self.i=nn.Linear(ATOM_DIM+BOND_DIM,h); self.u=nn.Linear(h,h,bias=False); self.o=nn.Linear(ATOM_DIM+h,h); self.steps=c.message_steps; self.d=nn.Dropout(c.dropout)
    def forward(self,b):
        if not b.edge_index.shape[1]: return torch.relu(self.o(torch.cat([b.atom_x,b.atom_x.new_zeros(len(b.atom_x),self.o.in_features-ATOM_DIM)],1)))
        s,t=b.edge_index; h0=torch.relu(self.i(torch.cat([b.atom_x[s],b.bond_x],1))); h=h0
        for _ in range(self.steps-1): inc=h.new_zeros(len(b.atom_x),h.shape[1]).index_add_(0,t,h); h=self.d(torch.relu(h0+self.u(inc[s]-h[b.reverse_edge])))
        inc=h.new_zeros(len(b.atom_x),h.shape[1]).index_add_(0,t,h); return torch.relu(self.o(torch.cat([b.atom_x,inc],1)))

class ChemicalEncoder(nn.Module):
    def __init__(self,c=None):
        super().__init__(); self.c=c or ChemicalEncoderConfig(); c=self.c; h=c.hidden_dim
        if c.hidden_dim%c.num_heads: raise ValueError("hidden_dim must be divisible by num_heads")
        self.m=_Morgan(c) if c.mode!="graph" else None; self.g=_DMPNN(c) if c.mode!="morgan" else None
        self.mg=nn.MultiheadAttention(h,c.num_heads,c.dropout,batch_first=True) if c.mode=="cross_attention" else None
        self.gm=nn.MultiheadAttention(h,c.num_heads,c.dropout,batch_first=True) if c.mode=="cross_attention" else None
        self.out=nn.Sequential(nn.Linear(h*(2 if c.mode in ("concat","cross_attention") else 1),c.output_dim),nn.LayerNorm(c.output_dim))
    def forward(self,b):
        mt=mm=gt=gm=mp=gp=None
        if self.m: mt,mm=self.m(b.morgan); mp=_mean(mt,mm)
        if self.g: gt,gm=_pad(self.g(b),b.atom_batch,len(b.morgan)); gp=_mean(gt,gm)
        if self.c.mode=="morgan": z=mp
        elif self.c.mode=="graph": z=gp
        elif self.c.mode=="concat": z=torch.cat([mp,gp],1)
        else:
            mc,_=self.mg(mt,gt,gt,key_padding_mask=gm); gc,_=self.gm(gt,mt,mt,key_padding_mask=mm)
            z=torch.cat([_mean(mt+mc,mm),_mean(gt+gc,gm)],1)
        return self.out(z)
