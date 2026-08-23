import argparse, json, random
from pathlib import Path
import numpy as np
import torch
from torch.utils.data import DataLoader, random_split
from .chemistry import ChemicalEncoderConfig
from .data import PerturbationDataset, collate_perturbations, write_demo_dataset
from .model import FusionCPA, FusionCPAConfig

def seed_all(seed): random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
def move(batch,device): return {k:(v.to(device) if hasattr(v,"to") else v) for k,v in batch.items()}
def main():
    p=argparse.ArgumentParser(); p.add_argument("--csv"); p.add_argument("--output",default="runs/demo"); p.add_argument("--mode",choices=["morgan","graph","concat","cross_attention"],default="cross_attention"); p.add_argument("--epochs",type=int,default=3); p.add_argument("--batch-size",type=int,default=16); p.add_argument("--seed",type=int,default=42); p.add_argument("--device",default="cuda" if torch.cuda.is_available() else "cpu"); a=p.parse_args(); seed_all(a.seed)
    out=Path(a.output); out.mkdir(parents=True,exist_ok=True); csv_path=a.csv or write_demo_dataset(out/"toy_data")
    ds=PerturbationDataset(csv_path); nval=max(1,len(ds)//5); train,val=random_split(ds,[len(ds)-nval,nval],generator=torch.Generator().manual_seed(a.seed))
    collate=lambda x:collate_perturbations(x); tl=DataLoader(train,a.batch_size,shuffle=True,collate_fn=collate); vl=DataLoader(val,a.batch_size,collate_fn=collate)
    genes=len(ds[0]["control"]); cfg=FusionCPAConfig(genes,len(ds.cell_vocab),chemical=ChemicalEncoderConfig(mode=a.mode)); model=FusionCPA(cfg).to(a.device); opt=torch.optim.AdamW(model.parameters(),lr=1e-3,weight_decay=1e-5); best=float("inf")
    for epoch in range(a.epochs):
        model.train(); tr=[]
        for b in tl:
            b=move(b,a.device); o=model(b["control"],b["molecules"],b["dose"],b["cell"]); loss=model.loss(o,b["target"]); opt.zero_grad(); loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(),5); opt.step(); tr.append(loss.item())
        model.eval(); va=[]
        with torch.no_grad():
            for b in vl: b=move(b,a.device); va.append(model.loss(model(b["control"],b["molecules"],b["dose"],b["cell"]),b["target"]).item())
        score=float(np.mean(va)); print(json.dumps({"epoch":epoch+1,"train_nll":float(np.mean(tr)),"val_nll":score}))
        if score<best: best=score; torch.save({"model":model.state_dict(),"config":vars(a),"cell_vocab":ds.cell_vocab,"n_genes":genes},out/"best.pt")
    with open(out/"metrics.json","w") as f: json.dump({"best_val_nll":best,"mode":a.mode,"seed":a.seed},f,indent=2)
if __name__=="__main__": main()

