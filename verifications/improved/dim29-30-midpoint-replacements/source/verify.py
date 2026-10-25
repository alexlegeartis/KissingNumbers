"""Portable exact verifier: NumPy, explicit coordinates, no search data.

Coordinates may have positive integer metric weights: physical coordinate j
is sqrt(weight[j])*x[j]. All head coordinates have weight one. New/new dot
products are integers with absolute partial sums <=32. BLAS float64 therefore
evaluates these integer operations exactly; --integer uses integer matmul.
No tolerance is used. Pure D5/E6/E7 root configurations are checked separately.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
import argparse,hashlib,itertools,json,time
from pathlib import Path
import numpy as np


def leech_shell():
    generator=np.array([1,0,1,0,1,1,1,0,0,0,1,1],dtype=np.uint8)
    matrix=np.zeros((12,23),dtype=np.uint8)
    for i in range(12):matrix[i,i:i+12]=generator
    bits=((np.arange(4096)[:,None]>>np.arange(12))&1).astype(np.uint8)
    golay=np.zeros((4096,24),dtype=np.uint8)
    golay[:,:23]=(bits@matrix)%2;golay[:,23]=golay[:,:23].sum(1)%2
    assert dict(zip(*np.unique(golay.sum(1),return_counts=True)))=={0:1,8:759,12:2576,16:759,24:1}
    rows=[]
    for i,j in itertools.combinations(range(24),2):
        for a,b in itertools.product([4,-4],repeat=2):
            row=np.zeros(24,dtype=np.int8);row[i]=a;row[j]=b;rows.append(row)
    signs=np.array([[(-2 if (k>>i)&1 else 2) for i in range(8)] for k in range(256) if bin(k).count('1')%2==0],dtype=np.int8)
    for word in golay[golay.sum(1)==8]:
        block=np.zeros((128,24),dtype=np.int8);block[:,np.flatnonzero(word)]=signs;rows.extend(block)
    sign_words=np.where(golay,-1,1).astype(np.int8)
    for j in range(24):
        block=sign_words.copy();block[:,j]*=-3;rows.extend(block)
    m=np.array(rows,dtype=np.int8)
    assert m.shape==(196560,24)
    assert np.all(np.sum(m.astype(np.int32)**2,axis=1)==32)
    assert len({row.tobytes() for row in m})==len(m)
    return m


def pure_roots(k):
    rows=[];base=k if k==5 else k-1
    for i,j in itertools.combinations(range(base),2):
        for a,b in itertools.product([-4,4],repeat=2):
            row=np.zeros(k,dtype=np.int16);row[i]=a;row[j]=b;rows.append(row)
    weight=np.ones(k,dtype=np.int16)
    if k==6:
        weight[-1]=3
        for signs in itertools.product([-1,1],repeat=5):
            rows.append(np.r_[2*np.array(signs),2*np.prod(signs)].astype(np.int16))
    elif k==7:
        weight[-1]=2
        for sign in [-1,1]:
            row=np.zeros(k,dtype=np.int16);row[-1]=4*sign;rows.append(row)
        for signs in itertools.product([-1,1],repeat=6):
            if np.prod(signs)!=1:continue
            for sign in [-1,1]:rows.append(np.r_[2*np.array(signs),2*sign].astype(np.int16))
    roots=np.array(rows,dtype=np.int16);g=(roots*weight)@roots.T
    assert np.all(np.diag(g)==32);np.fill_diagonal(g,-32);assert g.max()<=16
    assert len(roots)=={5:40,6:72,7:126}[k]
    return len(roots)


def verify(path,integer=False):
    start=time.time();data=np.load(path,allow_pickle=False)
    new=data['new_points'].astype(np.int16);removed=data['removed_leech'].astype(np.int16)
    endpoints=data['midpoint_owner_indices'].astype(np.int64);d=new.shape[1]
    weight=data['metric_weights'].astype(np.int16) if 'metric_weights' in data else np.ones(d,dtype=np.int16)
    assert d in [29,30,31] and weight.shape==(d,) and np.all(weight>0) and np.all(weight[:24]==1)
    m=leech_shell();shell={row.tobytes() for row in m}
    assert len(np.unique(removed,axis=0))==len(removed)
    assert np.all(abs(removed)<=4)
    assert all(row.astype(np.int8).tobytes() in shell for row in removed)
    assert endpoints.shape==(len(new),2) and np.all(endpoints>=0) and np.all(endpoints<len(removed))
    assert np.all(endpoints[:,0]!=endpoints[:,1])
    assert len(np.unique(new,axis=0))==len(new)
    assert np.all(np.sum(new.astype(np.int32)**2*weight,axis=1)==32)
    assert np.all(np.sum(new[:,24:].astype(np.int32)**2*weight[24:],axis=1)==8)
    assert np.all(2*new[:,:24]==removed[endpoints].sum(axis=1))
    # Each midpoint clears all retained Leech points, since both endpoints
    # are deleted members of the standard Leech shell.
    x=new if integer else new.astype(np.float64);worst=-100;contacts=0
    for at in range(0,len(x),512):
        block=x[at:at+512];dot=(block*weight)@x.T
        # Weighted Cauchy-Schwarz bounds every absolute partial sum by 32.
        # Every product and sum is an exactly represented small integer.
        assert np.all(dot==np.rint(dot))
        dot[np.arange(len(block)),at+np.arange(len(block))]=-100
        assert dot.max()<=16,(at,float(dot.max()))
        worst=max(worst,int(dot.max()));contacts+=int(np.sum(dot==16))
    pure=pure_roots(d-24)
    # Physical tail norm sqrt(8) and pure-root norm sqrt(32) imply product
    # <=16, independently of the relative coordinate fields/orientation.
    total=196560-len(removed)+len(new)+pure
    result=dict(status='PASS EXACT',dimension=d,points=total,removed=len(removed),inserted=len(new),pure=pure,
                maximum_inserted_product=worst,inserted_contacts=contacts//2,
                certificate_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),
                arithmetic='integer matmul' if integer else 'exact small integer operations in float64 BLAS; no tolerance',
                seconds=time.time()-start)
    print(json.dumps(result,indent=2));return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('certificate',type=Path);ap.add_argument('--integer',action='store_true')
    ap.add_argument('--output',type=Path);args=ap.parse_args();result=verify(args.certificate,args.integer)
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
