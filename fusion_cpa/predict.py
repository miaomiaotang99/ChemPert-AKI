import argparse, csv
from pathlib import Path
import numpy as np
import torch
from .chemistry import ChemicalEncoderConfig, collate_molecules
from .model import FusionCPA, FusionCPAConfig

def main():
    p=argparse.ArgumentParser(); p.add_argument("--checkpoint",required=True); p.add_argument("--control",required=True); p.add_argument("--smiles",required=True); p.add_argument("--dose",type=float,required=True); p.add_argument("--cell-type",required=True); p.add_argument("--output",required=True); p.add_argument("--device",default="cpu"); a=p.parse_args()
    ck=torch.load(a.checkpoint,map_location=a.device,weights_only=False); ca=ck["config"]; mode=ca["mode"]
    cfg=FusionCPAConfig(ck["n_genes"],len(ck["cell_vocab"]),chemical=ChemicalEncoderConfig(mode=mode)); model=FusionCPA(cfg).to(a.device); model.load_state_dict(ck["model"]); model.eval()
    if a.cell_type not in ck["cell_vocab"]: raise ValueError(f"Unknown cell type {a.cell_type!r}; choices: {sorted(ck['cell_vocab'])}")
    control=torch.from_numpy(np.load(a.control).astype("float32")).reshape(1,-1).to(a.device); mol=collate_molecules([a.smiles]).to(a.device)
    with torch.no_grad(): out=model(control,mol,torch.tensor([[a.dose]],device=a.device),torch.tensor([ck["cell_vocab"][a.cell_type]],device=a.device))
    np.savez_compressed(a.output,mean=out["mean"].cpu().numpy()[0],variance=out["variance"].cpu().numpy()[0])
    print(f"Saved prediction to {a.output}")
if __name__=="__main__": main()

