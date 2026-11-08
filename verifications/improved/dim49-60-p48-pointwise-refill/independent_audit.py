"""Independent read-only replay of antipodal pointwise-refill certificates.
Usage: python independent_audit.py certificate.npz [pole_count]
No optimizer IDs, cached conflicts, counts or receipts are trusted.
"""
import argparse,hashlib,io,json,sys,time
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits


def integers(a):
    assert a.dtype.kind in 'iu'
    assert a.size==0 or (int(a.min())>=-2**31 and int(a.max())<2**31)
    return a.astype(np.int64)

def maxabs(a): return max(abs(int(a.min())),abs(int(a.max()))) if a.size else 0

def product(a,b):
    # Inputs, products and every partial sum are exact integers in binary64.
    assert a.ndim==b.ndim==2 and a.shape[1]==b.shape[0]
    assert a.shape[1]*maxabs(a)*maxabs(b)<2**53
    return a.astype(np.float64)@b.astype(np.float64)

def canon(x):
    assert np.all(np.any(x!=0,axis=1))
    j=np.argmax(x!=0,axis=1)
    return x*np.sign(x[np.arange(len(x)),j])[:,None]

def normcheck(x,G):
    w=product(x,G)
    assert x.shape[1]*maxabs(w)*maxabs(x)<2**53
    assert np.all(np.sum(w*x,axis=1)==6)

def main(path,poles,source):
    start=time.monotonic();certificate_bytes=Path(path).read_bytes();a=np.load(io.BytesIO(certificate_bytes),allow_pickle=False)
    source=Path(source)
    family=source/'family_refilled.npz'
    gram=source/'p48p_gram.npy'
    if not family.exists():family=source/'october/family_refilled.npz'
    if not gram.exists():gram=source/'data/p48p_gram.npy'
    gram_bytes=gram.read_bytes();family_bytes=family.read_bytes()
    assert hashlib.sha256(gram_bytes).hexdigest()=='01a5d107e4424092e5300e201592aaff62804ee27eeb1d23b6e516249897d3e1', 'Published P48 Gram identity mismatch'
    z=np.load(io.BytesIO(family_bytes),allow_pickle=False)
    B=integers(z['base']);G=integers(np.load(io.BytesIO(gram_bytes),allow_pickle=False))
    assert B.shape==(7783,48) and G.shape==(48,48)
    assert np.array_equal(G,G.T) and np.array_equal(G,z['gram'])
    if 'base' in a:assert np.array_equal(a['base'],B)
    if 'gram' in a:assert np.array_equal(a['gram'],G)
    normcheck(B,G);BG=product(B,G)
    basepairs=0
    for s in range(0,len(B),256):
        H=product(BG[s:s+256],B.T)
        for j,row in enumerate(H):
            i=s+j
            assert np.all(abs(row[:i])<=2) and np.all(abs(row[i+1:])<=2)
            basepairs+=i
    T=integers(a['transforms']);m=len(T)
    assert T.shape==(m,48,48) and m>0
    bits=a['retained_bits'];assert bits.dtype==np.uint8 and bits.shape==(m,(len(B)+7)//8)
    bits=np.unpackbits(bits,axis=1);assert not bits[:,len(B):].any();bits=bits[:,:len(B)].astype(bool)
    E=integers(a['extra_points']);O=integers(a['extra_offsets'])
    assert E.ndim==2 and E.shape[1]==48 and O.shape==(m+1,)
    assert O[0]==0 and O[-1]==len(E) and np.all(np.diff(O)>=0)
    normcheck(E,G)
    allheads=[];sizes=[];extrapairs=0
    for i,t in enumerate(T):
        assert np.array_equal(product(product(t,G),t.T),G)
        # Materialize final physical coordinates before pair checks; do not reuse
        # optimizer's base-frame compatibility calculation.
        b=product(B[bits[i]],t);e=product(E[O[i]:O[i+1]],t)
        normcheck(b,G);normcheck(e,G)
        eg=product(e,G)
        if len(e) and len(b):assert np.all(abs(product(eg,b.T))<=2)
        extrapairs+=len(e)*len(b)
        if len(e):
            H=product(eg,e.T)
            assert np.all(np.diag(H)==6)
            np.fill_diagonal(H,0);assert np.all(abs(H)<=2)
            extrapairs+=len(e)*(len(e)-1)//2
        x=canon(np.concatenate([b,e]));assert maxabs(x)<2**15
        assert np.all(x==np.rint(x))
        allheads.append(x.astype(np.int16));sizes.append(len(x))
    X=np.concatenate(allheads)
    keys=np.ascontiguousarray(X).view(np.dtype((np.void,96))).ravel()
    assert len(np.unique(keys))==len(keys),'Duplicate antipodal line ownership'
    if 'sizes' in a:assert np.array_equal(a['sizes'],sizes)
    result=dict(status='PASS',certificate_sha256=hashlib.sha256(certificate_bytes).hexdigest(),classes=m,selected_lines=len(X),class_sizes=sizes,extra_lines=len(E),base_pairs=basepairs,extra_pairs=extrapairs,arithmetic='independent binary64 integer products with strict <2^53 absolute partial-sum bounds',disjointness='complete canonical 48-coordinate int16 records, no hashes',mother_shell_scope='Published integral P48 Gram lattice, minimum norm six and population 52416000 inherited; shell not enumerated',tail_pole_scope='Not checked by this head-only replay',elapsed_seconds=round(time.monotonic()-start,3))
    result['source_sha256']={'family_refilled.npz':hashlib.sha256(family_bytes).hexdigest(),'p48p_gram.npy':hashlib.sha256(gram_bytes).hexdigest()}
    assert Path(path).read_bytes()==certificate_bytes, 'Certificate changed during replay'
    assert gram.read_bytes()==gram_bytes and family.read_bytes()==family_bytes, 'Source changed during replay'
    if poles is not None:result['conditional_kissing_count']=52416000+2*len(X)+poles
    return result
if __name__=='__main__':
    if not __debug__:raise SystemExit('Run without -O')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate',type=Path)
    parser.add_argument('poles',type=int,nargs='?')
    parser.add_argument('--source',type=Path,default=Path(__file__).resolve().parent/'source',help='Directory containing family_refilled.npz and p48p_gram.npy; original cap-package root also accepted.')
    args=parser.parse_args()
    with threadpool_limits(limits=2):result=main(args.certificate,args.poles,args.source)
    print(json.dumps(result,indent=2))
