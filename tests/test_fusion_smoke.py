import pytest, torch
from fusion_cpa import ChemicalEncoderConfig, FusionCPA, FusionCPAConfig, collate_molecules

@pytest.mark.parametrize("mode",["morgan","graph","concat","cross_attention"])
def test_end_to_end(mode):
    c=ChemicalEncoderConfig(mode=mode,hidden_dim=32,output_dim=16,num_heads=4,dropout=0)
    model=FusionCPA(FusionCPAConfig(n_genes=20,n_cell_types=2,latent_dim=16,hidden_dim=32,chemical=c))
    out=model(torch.randn(3,20),collate_molecules(["CCO","c1ccccc1","[Na+]"]),torch.rand(3,1),torch.tensor([0,1,0]))
    assert out["mean"].shape==out["variance"].shape==(3,20); model.loss(out,torch.randn(3,20)).backward()

