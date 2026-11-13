"""Cancel unchanged rows and independently replay the30D section exchange."""
from pathlib import Path
import hashlib,json,time
import numpy as np


def main():
    t0=time.monotonic();here=Path(__file__).resolve().parent;out=here/'ls15_section30';root=here.parents[2]
    src=root.parent/'PackingStar/25D-31D/new_bounds_20261111/D30_247992_coordinates.npz'
    assert hashlib.sha256(src.read_bytes()).hexdigest()=='5e46243e5f3cf4790d56abc4008c9cebeb45bdde2fd3080e44a0546d59cccc5a'
    z=np.load(src);X=2*z['points'].astype(np.int16);assert int(z['coordinate_denominator'])==1
    M=np.load(here.parent/'algebraic/ls15_exact/module.npz')['points'].astype(np.int64)
    section=M[M[:,14]==0][:,:14];r=np.load(out/'best_module.npz');R=r['transformation_numerator'].astype(np.int64);D=int(r['transformation_denominator'])
    assert D==4 and np.array_equal(R@R.T,16*np.eye(14,dtype=np.int64))
    num=section@R.T;assert np.all(num%4==0);Y=num//4;assert np.array_equal(Y,r['points']) and len(Y)==1356
    assert np.all(np.sum(Y*Y,axis=1)==32)
    G=Y@Y.T;np.fill_diagonal(G,-32);assert G.max()<=16
    embed=np.r_[[0,2,4,5,6,10,11,23],np.arange(24,30)]
    profiles,inverse=np.unique(X[:,embed],axis=0,return_inverse=True)
    blocked=np.any((2*Y).astype(np.float32)@profiles.astype(np.float32).T>64,axis=0)
    removed=X[blocked[inverse]]
    added=np.zeros((len(Y),30),np.int16);added[:,embed]=2*Y
    assert len(removed)==1264 and len(added)==1356
    rd={v.tobytes():i for i,v in enumerate(removed)};ad={v.tobytes():i for i,v in enumerate(added)};common=rd.keys()&ad.keys()
    removed=removed[[i for k,i in rd.items() if k not in common]];added=added[[i for k,i in ad.items() if k not in common]]
    np.savez_compressed(out/'exchange_minimal.npz',removed=removed,added=added)
    xd={v.tobytes():i for i,v in enumerate(X)};assert all(v.tobytes() in xd for v in removed)
    keep=np.ones(len(X),bool);keep[[xd[v.tobytes()] for v in removed]]=False
    F=np.vstack([X[keep],added]);expected=np.load(out/'coordinates30_248084.npz')['points'].astype(np.int16)
    assert len(F)==248084 and np.array_equal(np.unique(F,axis=0),np.unique(expected,axis=0))
    assert np.all(np.sum(F.astype(np.int64)**2,axis=1)==128)
    retained=X[keep].astype(np.int64);maximum=-128
    for row in added.astype(np.int64):
        dots=np.zeros(len(retained),np.int64)
        for j in np.flatnonzero(row):dots+=int(row[j])*retained[:,j]
        assert dots.max(initial=-128)<=64;maximum=max(maximum,int(dots.max(initial=-128)))
    gg=added.astype(np.int64)@added.astype(np.int64).T;np.fill_diagonal(gg,-128);assert gg.max(initial=-128)<=64
    heads=Y[np.sum(Y[:,8:]**2,axis=1)==16];tails,cnt=np.unique(heads[:,8:],axis=0,return_counts=True);vals,counts=np.unique(cnt,return_counts=True)
    report=dict(status='PASS exact transform, minimal exchange and independent direct int64 changed pairs',baseline=len(X),points=len(F),gain=len(F)-len(X),
                canonical_section_omitted_tail=6,transformation_denominator=D,module_points=len(Y),module_heads=len(heads),
                tail_fiber_occupancies=dict(zip(map(str,vals),map(int,counts))),
                whole_removed=1264,whole_added=1356,common_cancelled=len(common),removed=len(removed),added=len(added),
                new_retained_pairs=len(added)*len(retained),new_new_pairs=len(added)*(len(added)-1)//2,
                maximum_new_retained_product=maximum,squared_norm=128,inner_product_limit=64,
                source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),coordinates_sha256=hashlib.sha256((out/'coordinates30_248084.npz').read_bytes()).hexdigest(),
                minimal_exchange_sha256=hashlib.sha256((out/'exchange_minimal.npz').read_bytes()).hexdigest(),
                seconds=time.monotonic()-t0,scope='Changed geometry independently checked. Unchanged source pairs require the separate complete final verifier.')
    (out/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
