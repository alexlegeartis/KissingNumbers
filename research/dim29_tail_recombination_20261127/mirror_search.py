#!/usr/bin/env python3
"""Exact stock-preserving mixing of two isometric Takhanov v0.2 29D copies."""
import hashlib, itertools, json, pathlib, sys
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching

SOURCE_BLOB = "8e9dd3f474f8247214fe6cbe9fd00c15c1ddb07a"
here = pathlib.Path(__file__).resolve().parent
z = (here / "d29_certificate.npz").read_bytes()
assert hashlib.sha1(f"blob {len(z)}\0".encode()+z).hexdigest()==SOURCE_BLOB
q = np.load(here / "d29_certificate.npz", allow_pickle=False)
stock = q["added_integer"].astype("<i2")
assert stock.shape==(67868,29)
assert np.all(np.einsum("ij,ij->i",stock.astype(np.int32),stock.astype(np.int32))==128)
assert np.all(np.abs(stock) <= 11)
def key(a): return a.copy(order="C").view("V58").reshape(-1)
stock_keys=key(stock)
assert len(np.unique(stock_keys))==len(stock)
assert (len(stock)-len(q["deleted_leech_indices"])+196560)==232230

def run(label, perm, signs):
    other=stock.copy()
    other[:,24:]=stock[:,24+np.asarray(perm)]*np.asarray(signs)
    shared,i,j=np.intersect1d(stock_keys,key(other),assume_unique=True,return_indices=True)
    oa=np.ones(len(stock),dtype=bool);ob=oa.copy()
    oa[i]=False; ob[j]=False
    left=other[ob]; right=stock[oa]
    n=len(left)
    if n == 0:
        return dict(label=label, asymmetric=0, gain=0)
    if n>48000:
        return dict(label=label, asymmetric=n, skipped=True)
    # Integer bounded partial sums guarantee these binary32 products exact.
    lhs=left.astype(np.float32); rhs=right.astype(np.float32).T.copy()
    rows=[];cols=[]
    for a in range(0,n,256):
        ii,jj=np.nonzero(lhs[a:a+256]@rhs > 64)
        rows.append(ii+a);cols.append(jj)
    rows=np.concatenate(rows); cols=np.concatenate(cols)
    G=csr_matrix((np.ones(len(rows),dtype=np.int8),(rows,cols)),shape=(n,n))
    match_left=maximum_bipartite_matching(G,perm_type="column")
    assert len(match_left)==n
    gain=int(np.count_nonzero(match_left<0))
    out=dict(label=label,asymmetric=n,edges=len(rows),gain=gain)
    if not gain:
        assert len(np.unique(match_left))==n
        assert np.all(np.einsum("ij,ij->i",left.astype(np.int32),right[match_left].astype(np.int32))>64)
        np.save(here/(label+"_matching.npy"),match_left.astype("<i4"))
    if gain:
        # Alternating search returns a largest Hall-deficient donor set.
        match_right=np.full(n,-1)
        for a,b in enumerate(match_left):
            if b>=0: match_right[b]=a
        J=set(np.flatnonzero(match_left<0)); K=set(); work=list(J)
        while work:
            a=work.pop()
            for b in G.indices[G.indptr[a]:G.indptr[a+1]]:
                if b in K: continue
                K.add(int(b)); c=int(match_right[b])
                if c>=0 and c not in J: J.add(c); work.append(c)
        u=np.array(sorted(J),dtype=np.int32); v=np.array(sorted(K),dtype=np.int32)
        assert len(u)-len(v)==gain
        d=left[u]; delete=np.flatnonzero(oa)[v]
        keep=np.ones(len(stock),dtype=bool);keep[delete]=False
        for p in range(0,len(d),128):
            assert np.max(d[p:p+128].astype(np.int32)@stock[keep].astype(np.int32).T)<=64
        np.savez_compressed(here / (label+"_gain.npz"),added=d,deleted_indices=delete,
                            tail_permutation=perm,tail_signs=signs)
        out.update(added=len(d),removed=len(delete))
    return out

p=tuple(range(5))
trials=[("negative_tail",p,(-1,)*5)]
if len(sys.argv)>1:
    trials.extend([("flip_0",p,(-1,1,1,1,1)),("swap_1_2",(0,2,1,3,4),(1,)*5),
                   ("quarter_0_4",(4,1,2,3,0),(-1,1,1,1,1))])
for trial in trials:
    print(json.dumps(run(*trial)),flush=True)
