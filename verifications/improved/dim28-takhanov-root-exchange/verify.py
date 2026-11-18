"""Standalone exact verification of the frozen 28D continuation.

No search program, external checkout, solver or lattice premise is used.
Every final unordered pair is covered. A separate int64 pass verifies the
source exchange. Binary32 is used only for bounded integer products, with
an explicit exact-representability bound and no numerical tolerance.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import time
import numpy as np

ROOT=Path(__file__).resolve().parent


def need(condition,message):
    if not condition: raise ValueError(message)


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(X):
    return X[np.lexsort(X.T[::-1])].astype('<i2')


def row_hash(X):
    return hashlib.sha256(canonical(X).tobytes()).hexdigest()


def integer_geometry(points,block=384,direct=False):
    need(points.ndim==2 and len(points)>0,'coordinate shape')
    need(np.issubdtype(points.dtype,np.integer),'coordinates must be integers')
    X=points.astype(np.int64)
    largest=int(np.max(np.abs(X)))
    bound=X.shape[1]*largest*largest
    need(bound<2**63,'int64 norm overflow')
    need(np.all(np.einsum('ij,ij->i',X,X)==128),'squared norm must equal128')
    need(len(np.unique(X,axis=0))==len(X),'duplicate vector')
    # A repeated canonical representative is a genuine antipodal pair because
    # original duplicates were already rejected. Singles keep their actual sign.
    first=np.argmax(X!=0,axis=1)
    sign=np.sign(X[np.arange(len(X)),first])
    need(np.all(sign!=0),'zero vector')
    C=X*sign[:,None]
    _,where,multiplicity=np.unique(C,axis=0,return_index=True,return_counts=True)
    need(np.all((multiplicity==1)|(multiplicity==2)),'bad antipodal multiplicity')
    P=C[where[multiplicity==2]]
    U=X[where[multiplicity==1]]
    dtype=np.int64 if direct else np.float32
    need(bound<(2**63 if direct else 2**24),'matrix-product exactness bound')
    peak=-128
    checked_kernels=0
    for A,B,absolute,triangle in [(P,P,True,True),(P,U,True,False),(U,U,False,True)]:
        if not len(A) or not len(B):continue
        AA,BB=A.astype(dtype),B.astype(dtype)
        for at in range(0,len(A),block):
            stop=min(len(A),at+block)
            products=AA[at:stop]@BB[at if triangle else 0:].T
            if absolute:np.abs(products,out=products)
            if triangle:
                diagonal=np.arange(stop-at)
                products[diagonal,diagonal]=-128
            value=int(products.max())
            need(value<=64,'distinct-point inner product exceeds64')
            peak=max(peak,value)
            checked_kernels+=1
    p,u=len(P),len(U)
    covered=4*math.comb(p,2)+p+2*p*u+math.comb(u,2)
    need(covered==math.comb(len(X),2),'unordered-pair coverage')
    return dict(count=len(X),dimension=X.shape[1],antipodal_pairs=p,unpaired_points=u,
                unordered_pairs=covered,maximum_inner_product=peak,
                squared_norm=128,inner_product_limit=64,
                kernel='direct int64' if direct else 'binary32 exact bounded integers',
                absolute_partial_sum_bound=bound,matrix_blocks=checked_kernels)


def exchange_check(source,final,removed,added):
    need(source.shape[1]==final.shape[1]==removed.shape[1]==added.shape[1],'exchange dimension')
    need(np.issubdtype(removed.dtype,np.integer)and np.issubdtype(added.dtype,np.integer),'exchange integer dtype')
    src={r.tobytes():i for i,r in enumerate(source.astype('<i2'))}
    need(len(src)==len(source),'duplicate source point')
    minus=[r.tobytes()for r in removed.astype('<i2')]
    need(len(set(minus))==len(minus)and all(r in src for r in minus),'removed points must be distinct source members')
    keep=np.ones(len(source),bool)
    keep[[src[r]for r in minus]]=False
    retained=source[keep]
    rebuilt=np.r_[retained,added]
    need(np.array_equal(canonical(rebuilt),canonical(final)),'source exchange reconstruction')
    need(len(final)==len(source)-len(removed)+len(added),'exchange count')
    # An arithmetically independent changed-pair proof uses signed int64.
    A=added.astype(np.int64);R=retained.astype(np.int64)
    need(np.all(np.einsum('ij,ij->i',A,A)==128),'added-point norm')
    peak=-128
    for at in range(0,len(A),16):
        products=A[at:at+16]@R.T
        if products.size:peak=max(peak,int(products.max()))
        need(not products.size or int(products.max())<=64,'added/retained conflict')
    for i in range(len(A)):
        products=A[i+1:]@A[i]
        if products.size:peak=max(peak,int(products.max()))
        need(not products.size or int(products.max())<=64,'added/added conflict')
    return dict(source_count=len(source),retained=len(retained),removed=len(removed),added=len(added),
                gain=len(added)-len(removed),changed_unordered_pairs=len(A)*len(R)+math.comb(len(A),2),
                maximum_changed_inner_product=peak,kernel='direct signed int64')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--report',type=Path);ap.add_argument('--direct',action='store_true');a=ap.parse_args();start=time.monotonic()
    manifest=json.loads((ROOT/'manifest.json').read_text())
    for name,digest in manifest['files'].items():
        need(file_hash(ROOT/name)==digest,'file checksum: '+name)
    data=ROOT/'data'
    with np.load(ROOT/manifest['coordinate_file'],allow_pickle=False)as z:
        X=z['points'];need(int(z['squared_norm'])==128 and int(z['inner_product_limit'])==64,'coordinate metadata')
    need(X.shape==(manifest['count'],manifest['dimension']),'frozen count and dimension')
    need(row_hash(X)==manifest['canonical_coordinate_sha256'],'canonical coordinate identity')
    with np.load(ROOT/manifest['source_file'],allow_pickle=False)as z:S=z['points']
    with np.load(data/'exchange.npz',allow_pickle=False)as z:removed,added=z['removed'],z['added']
    delta=exchange_check(S,X,removed,added)
    print('EXACT_EXCHANGE',json.dumps(delta),flush=True)
    geometry=integer_geometry(X,direct=a.direct)
    report=dict(verified=True,complete_final_geometry=True,independent_changed_pair_check=True,
                geometry=geometry,exchange=delta,coordinate_file=manifest['coordinate_file'],
                canonical_coordinate_sha256=row_hash(X),source=manifest['source'],
                seconds=time.monotonic()-start)
    if a.report is not None:
        a.report.write_text(json.dumps(report,indent=2)+'\n')
    print('VERIFIED',json.dumps(report),flush=True)


if __name__=='__main__':main()
