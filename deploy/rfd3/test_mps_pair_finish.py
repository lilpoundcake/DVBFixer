import pytest
import torch
from mps_pair_finish import finish
from rfd3.model.layers.block_utils import pairwise_mean_pool


@pytest.mark.parametrize('device',['cpu','mps'])
def test_fused_initializer_vs_native(device):
    idx = torch.tensor([0,1,1,1,2,2,4,4,4,4,5],device=device)
    m = torch.nn.Module().eval()
    for name in ('process_single_l','process_single_m','process_z','process_pll'):
        setattr(m,name,torch.nn.Linear(16,16).to(device))
    m.pair_mlp=torch.nn.Sequential(torch.nn.ReLU(),torch.nn.Linear(16,16)).to(device)
    pair=torch.randn(11,11,16,device=device)
    single=torch.randn(11,16,device=device)
    tokens=torch.randn(6,6,16,device=device)
    with torch.no_grad():
        expected=pair+(m.process_single_l(single)[:,None]+m.process_single_m(single)[None,:])
        expected=expected+m.process_z(tokens)[idx][:,idx]
        expected=expected+m.pair_mlp(expected)
        pool=pairwise_mean_pool(m.process_pll(expected)[None],idx,6,pair.dtype)[0]
        actual=finish(m,pair,single,tokens,idx,rows=3)
    torch.testing.assert_close(expected,pair)
    torch.testing.assert_close(actual,pool)


@pytest.mark.parametrize('device',['cpu','mps'])
def test_local_bias_reordering(device):
    from mps_pair_finish import gather_pair, sparse_rows
    from rfd3.model.layers.attention import sparse_pairbias_attention
    q,k,v=torch.randn(3,2,11,16,device=device).unbind(0)
    idx=torch.randint(11,(2,11,4),device=device)
    proj=torch.nn.Linear(16,4).to(device)
    for shape in [(11,11,16),(2,11,11,16)]:
        pair=torch.randn(shape,device=device)
        with torch.no_grad():
            baseline=sparse_pairbias_attention(q,k,v,proj(pair),idx,4)
            actual=sparse_pairbias_attention(q,k,v,proj(gather_pair(pair,idx)),idx,4,gather_bias=False)
        torch.testing.assert_close(actual,baseline)
        with torch.no_grad():
            gate=torch.randn_like(q)
            tiled=sparse_rows(q,k,v,proj(gather_pair(pair,idx)),idx,4,G=gate,rows=3)
            native=sparse_pairbias_attention(q,k,v,proj(pair),idx,4,G=gate)
        torch.testing.assert_close(tiled,native)
