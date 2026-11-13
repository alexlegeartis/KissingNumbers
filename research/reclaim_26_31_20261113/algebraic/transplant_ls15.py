"""Place a verified 2564-point E8-equator module in the complete 31D source."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np

NORMAL=np.array([0,2,4,5,6,10,11,23])
EMBED=np.r_[NORMAL,np.arange(24,31)]


def transplant(source,module,out):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=True)
    z=np.load(source,allow_pickle=False)
    assert int(z['squared_norm'])==128 and np.all(z['points']%2==0)
    x=z['points'].astype(np.int16)//2
    m=np.load(module,allow_pickle=False)['points'].astype(np.int16)
    assert len(m)==2564 and np.all(np.sum(m.astype(np.int32)**2,axis=1)==32)
    gm=m.astype(np.int32)@m.T.astype(np.int32);np.fill_diagonal(gm,-32)
    assert gm.max()<=16
    outside=[j for j in range(31) if j not in EMBED]
    oldmodule=~np.any(x[:,outside],axis=1)
    oldeq=x[oldmodule&~np.any(x[:,24:],axis=1)][:,NORMAL]
    eq=~np.any(m[:,8:],axis=1)
    target={r.tobytes() for r in oldeq}
    flip=False
    if {r[:8].tobytes() for r in m[eq]}!=target:
        m[:,7]*=-1;flip=True
    assert {r[:8].tobytes() for r in m[eq]}==target
    y=np.zeros((len(m),31),dtype=np.int16);y[:,EMBED]=m
    proj,inv,counts=np.unique(x[:,EMBED],axis=0,return_inverse=True,return_counts=True)
    # All products and every partial sum are exact bounded integers in binary32.
    assert 31*16*16 < 2**24
    blocked=np.any(m.astype(np.float32)@proj.astype(np.float32).T>16,axis=0)
    remove=blocked[inv]
    f=np.r_[x[~remove],y]
    assert len(f)==len(np.unique(f,axis=0))
    np.savez_compressed(out/f'coordinates31_{len(f)}.npz',points=2*f,squared_norm=128,
                        inner_product_limit=64,coordinate_denominator=2)
    np.savez_compressed(out/'exchange.npz',removed=2*x[remove],added=2*y)
    np.savez_compressed(out/'oriented_module.npz',points=m,squared_norm=32,inner_product_limit=16)
    # Independent direct int64 replay against every kept source point.
    checks=0;maximum=-32
    for row in y:
        values=np.zeros(np.sum(~remove),dtype=np.int64)
        for j in np.flatnonzero(row):values+=int(row[j])*x[~remove,j].astype(np.int64)
        assert values.max(initial=-32)<=16
        maximum=max(maximum,int(values.max(initial=-32)));checks+=len(values)
    report=dict(source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                module=str(module),module_sha256=hashlib.sha256(module.read_bytes()).hexdigest(),
                baseline=len(x),points=len(f),gain=len(f)-len(x),module_points=len(m),
                removed=int(np.sum(remove)),removed_module=int(np.sum(remove&oldmodule)),
                removed_external=int(np.sum(remove&~oldmodule)),old_module=int(np.sum(oldmodule)),
                flipped_core_coordinate_7=flip,maximum_new_old_product=maximum,
                independent_direct_int64_pairs=checks,
                new_new_unordered_pairs=len(m)*(len(m)-1)//2,
                changed_geometry='PASS',seconds=time.monotonic()-started,
                scope='Complete changed-pair proof. Unchanged source pairs inherit the pinned source until an independent full replay.')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('source',type=Path);ap.add_argument('module',type=Path);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    transplant(a.source,a.module,a.out)
