"""chemCPA-style latent intervention model with online chemical encoding."""
from dataclasses import dataclass, field
import torch
from torch import nn
from .chemistry import ChemicalEncoder, ChemicalEncoderConfig

@dataclass
class FusionCPAConfig:
    n_genes:int; n_cell_types:int; latent_dim:int=256; hidden_dim:int=512
    chemical:ChemicalEncoderConfig=field(default_factory=ChemicalEncoderConfig)

def _mlp(dims,dropout=.1):
    layers=[]
    for a,b in zip(dims[:-2],dims[1:-1]): layers += [nn.Linear(a,b),nn.LayerNorm(b),nn.GELU(),nn.Dropout(dropout)]
    return nn.Sequential(*layers,nn.Linear(dims[-2],dims[-1]))

class FusionCPA(nn.Module):
    """control expression + SMILES + dose + cell -> Gaussian perturbed expression."""
    def __init__(self,c):
        super().__init__(); self.config=c; c.chemical.output_dim=c.latent_dim
        self.chemical=ChemicalEncoder(c.chemical); self.cell=nn.Embedding(c.n_cell_types,c.latent_dim)
        self.encoder=_mlp([c.n_genes,c.hidden_dim,c.latent_dim]); self.doser=_mlp([1,64,c.latent_dim])
        self.decoder=_mlp([c.latent_dim,c.hidden_dim,2*c.n_genes])
    def forward(self,control,molecules,dose,cell):
        z0=self.encoder(control); zchem=self.chemical(molecules)
        gate=torch.sigmoid(self.doser(torch.log1p(dose.clamp_min(0))))
        z=z0+gate*zchem+self.cell(cell); stats=self.decoder(z); mean,raw_var=stats.chunk(2,1)
        return {"mean":mean,"variance":nn.functional.softplus(raw_var)+1e-4,"basal":z0,"chemical":zchem,"treated":z}
    @staticmethod
    def loss(output,target):
        return .5*(torch.log(output["variance"])+(target-output["mean"]).square()/output["variance"]).mean()

