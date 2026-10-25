"""Expand the two final certificates into complete exact coordinate arrays.

The output is (a + b*sqrt(3)), with squared norm 32. Divide by sqrt(32)
for unit vectors. This script supports the final 29D/30D certificates,
whose inserted-point metric weights are all one.
"""
import argparse,itertools
from pathlib import Path
import numpy as np
from verify import leech_shell


def materialize(certificate,output):
    data=np.load(certificate,allow_pickle=False)
    new=data['new_points'];removed=data['removed_leech'];dimension=new.shape[1]
    assert dimension in [29,30] and np.all(data['metric_weights']==1)
    deleted={r.tobytes() for r in removed.astype(np.int8)}
    shell=leech_shell();keep=np.array([r.tobytes() not in deleted for r in shell])
    old=np.pad(shell[keep],((0,0),(0,dimension-24)))
    pure_a=[];pure_b=[]
    for i,j in itertools.combinations(range(5),2):
        for s,t in itertools.product([-4,4],repeat=2):
            row=np.zeros(dimension,dtype=np.int8);row[24+i]=s;row[24+j]=t
            pure_a.append(row);pure_b.append(np.zeros(dimension,dtype=np.int8))
    if dimension==30:
        for signs in itertools.product([-1,1],repeat=5):
            a=np.zeros(dimension,dtype=np.int8);a[24:29]=2*np.array(signs)
            b=np.zeros(dimension,dtype=np.int8);b[29]=2*np.prod(signs)
            pure_a.append(a);pure_b.append(b)
    a=np.concatenate([old,new,np.array(pure_a)],axis=0)
    b=np.concatenate([np.zeros_like(old),np.zeros_like(new),np.array(pure_b)],axis=0)
    aa=a.astype(np.int32);bb=b.astype(np.int32)
    assert np.all((aa*aa+3*bb*bb).sum(1)==32) and np.all((aa*bb).sum(1)==0)
    assert len(np.unique(np.concatenate([a,b],axis=1),axis=0))==len(a)
    np.savez_compressed(output,a=a,b=b)
    print('Expanded',len(a),'exact points in dimension',dimension,'to',output)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('certificate',type=Path);ap.add_argument('output',type=Path)
    args=ap.parse_args();materialize(args.certificate,args.output)
