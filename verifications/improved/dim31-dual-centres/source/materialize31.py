"""Expand to exact (a+b*sqrt(2)+c*sqrt(3)+d*sqrt(6))/3 coordinates."""
import argparse
from pathlib import Path
import numpy as np
from verify31 import leech_shell


def materialize(path,output):
    z=np.load(path,allow_pickle=False);deleted={row.tobytes() for row in z['removed_leech'].astype(np.int8)}
    shell=leech_shell();bulk=shell[[row.tobytes() not in deleted for row in shell]]
    mid=z['midpoints'];centre=z['centres'];pure=z['pure_tails'];n=len(bulk)+len(mid)+len(centre)+len(pure)
    a=np.zeros((n,31),dtype=np.int8);b=np.zeros_like(a);c=np.zeros_like(a);d=np.zeros_like(a);at=0
    a[:len(bulk),:24]=3*bulk;at+=len(bulk)
    a[at:at+len(mid),:30]=3*mid[:,:30];b[at:at+len(mid),30]=3*mid[:,30];at+=len(mid)
    scale=tuple(map(int,z['centre_scale_squared'])) if 'centre_scale_squared' in z else (2,1)
    a[at:at+len(centre),:24]=3*centre[:,:24]
    if scale==(2,1):
        b[at:at+len(centre),24:30]=3*centre[:,24:30];a[at:at+len(centre),30]=6*centre[:,30]
    else:
        assert scale==(8,3)
        d[at:at+len(centre),24:30]=2*centre[:,24:30];c[at:at+len(centre),30]=4*centre[:,30]
    at+=len(centre);a[at:,24:30]=3*pure[:,:6];b[at:,30]=3*pure[:,6]
    aa,bb,cc,dd=[v.astype(np.int32) for v in [a,b,c,d]]
    assert np.all((aa*aa+2*bb*bb+3*cc*cc+6*dd*dd).sum(1)==288)
    assert not np.any((aa*bb+3*cc*dd).sum(1))
    assert not np.any((aa*cc+2*bb*dd).sum(1))
    assert not np.any((aa*dd+bb*cc).sum(1))
    assert len(np.unique(np.c_[a,b,c,d],axis=0))==n
    np.savez_compressed(output,a=a,b=b,c=c,d=d,denominator=np.array(3))
    print('Expanded',n,'exact 31D points to',output,flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('certificate',type=Path);ap.add_argument('output',type=Path)
    args=ap.parse_args();materialize(args.certificate,args.output)
