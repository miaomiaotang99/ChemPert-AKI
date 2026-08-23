"""CSV contract and deterministic splits for paired expression perturbations."""
import csv, json
from pathlib import Path
import numpy as np
import torch
from torch.utils.data import Dataset
from .chemistry import collate_molecules

class PerturbationDataset(Dataset):
    """CSV columns: smiles,dose,cell_type,control_path,perturbed_path.

    Expression files are .npy float vectors with an identical gene order.
    """
    def __init__(self,csv_path,cell_vocab=None):
        root=Path(csv_path).resolve().parent
        with open(csv_path,newline="",encoding="utf8") as f: self.rows=list(csv.DictReader(f))
        required={"smiles","dose","cell_type","control_path","perturbed_path"}
        if not self.rows or not required.issubset(self.rows[0]): raise ValueError(f"CSV requires {sorted(required)}")
        cells=sorted({r["cell_type"] for r in self.rows}); self.cell_vocab=cell_vocab or {x:i for i,x in enumerate(cells)}; self.root=root
    def __len__(self): return len(self.rows)
    def __getitem__(self,i):
        r=self.rows[i]
        return {"control":np.load(self.root/r["control_path"]).astype("float32"),"target":np.load(self.root/r["perturbed_path"]).astype("float32"),"smiles":r["smiles"],"dose":float(r["dose"]),"cell":self.cell_vocab[r["cell_type"]]}

def collate_perturbations(items,n_bits=2048,radius=2):
    return {"control":torch.from_numpy(np.stack([x["control"] for x in items])),"target":torch.from_numpy(np.stack([x["target"] for x in items])),"dose":torch.tensor([x["dose"] for x in items]).unsqueeze(1),"cell":torch.tensor([x["cell"] for x in items]),"molecules":collate_molecules([x["smiles"] for x in items],radius,n_bits)}

def write_demo_dataset(directory,n=64,genes=128,seed=42):
    """Self-contained nonlinear toy data for installation smoke tests only."""
    p=Path(directory); p.mkdir(parents=True,exist_ok=True); rng=np.random.default_rng(seed); smiles=["CCO","CC(=O)O","c1ccccc1","CCN"]
    rows=[]
    for i in range(n):
        c=rng.normal(size=genes).astype("float32"); dose=float(rng.uniform(.01,10)); j=i%4; y=(c+.15*np.log1p(dose)*(j+1)+rng.normal(0,.03,genes)).astype("float32")
        np.save(p/f"c{i}.npy",c); np.save(p/f"y{i}.npy",y); rows.append({"smiles":smiles[j],"dose":dose,"cell_type":f"cell_{i%2}","control_path":f"c{i}.npy","perturbed_path":f"y{i}.npy"})
    with open(p/"samples.csv","w",newline="",encoding="utf8") as f: w=csv.DictWriter(f,fieldnames=rows[0]); w.writeheader(); w.writerows(rows)
    with open(p/"metadata.json","w") as f: json.dump({"genes":genes,"seed":seed,"warning":"synthetic smoke-test data"},f)
    return p/"samples.csv"

