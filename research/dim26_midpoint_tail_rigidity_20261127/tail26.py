from collections import Counter
import hashlib,itertools,json,sys,pathlib
import numpy as np
from scipy.sparse import csr_matrix,bmat
from scipy.sparse.csgraph import connected_components
here=pathlib.Path(__file__).resolve().parent
b=(here/"common.py").read_bytes()
assert hashlib.sha1(f"blob {len(b)}\0".encode()+b).hexdigest()=="10d8321d8c0f831d474c3769adae7cc7caf29821"
from common import leech_shell
z=here/"d26_certificate.npz"
assert hashlib.sha256(z.read_bytes()).hexdigest()=="35fd8143cf0c7e6b0146719485e00c86952cbd5e404c7881c5f6b92fa25b21ec"
doc=np.load(z,allow_pickle=False)
A=doc["added_integer"].astype(np.int32)
M=doc["midpoint_rows"]; FR=doc["fraction_numerators"];FD=doc["fraction_denominators"]
R=A[np.setdiff1d(np.arange(len(A)),M)]
P,ind,ct=np.unique(A[M,:24],axis=0,return_index=True,return_counts=True)
assert P.shape==(4416,24) and (ct==2).all()
Q=np.array([[4,4],[4,-4]],np.int32)
pcol=np.argmax(abs(A[M,24:][ind]@Q.T),axis=1)
assert (np.einsum("ij,ij->i",P,P)==96).all()
old_pairs=set(map(tuple,A[M].tolist()))
assert all(tuple(np.r_[b,t]) in old_pairs for b,c in zip(P,pcol) for t in (Q[c],-Q[c]))
K=np.hstack((R[:,:24].astype(np.int64),R[:,24:].astype(np.int64)))
FB=np.vstack((K[:,:24],FR[:,:24]));FT=np.vstack((K[:,24:],FR[:,24:]))
Fden=np.r_[np.ones(len(R),np.int64),FD]

def frozen_constraints(U2): # input double body, output [n,2] true for feasible paired root line
    v=U2@FB.T
    return np.column_stack([np.max(v+2*np.abs(FT@t)-128*Fden,axis=1)<=0 for t in Q])

def bank_graph():
    X,Y=np.flatnonzero(pcol==0),np.flatnonzero(pcol==1)
    e=csr_matrix(P[X]@P[Y].T>32)
    c,m=connected_components(bmat([[None,e],[e.T,None]],format="csr"),directed=False)
    g=np.empty(len(P),np.int32);g[np.r_[X,Y]]=m
    fixed_ind=np.flatnonzero(~frozen_constraints(2*P)[np.arange(len(P)),1-pcol])
    pinned=set(map(int,g[fixed_ind]))
    return g,pinned,c

def candidates():
    shell=leech_shell().astype(np.int32);D=doc["deleted_leech_indices"];dead=shell[D]
    known=set(map(tuple,(2*P).tolist()))
    which={}
    for i in range(len(dead)):
        v=dead@dead[i]
        for j in np.flatnonzero(v[i+1:]==64)+i+1:
            b=tuple(map(int,dead[i]+dead[j]))
            if b not in known:which.setdefault(b,(int(D[i]),int(D[j])))
    return which

if __name__=="__main__":
    g,pinned,n=bank_graph()
    print(json.dumps(dict(components=n,pinned=len(pinned))),flush=True)
    C=candidates()
    print(json.dumps(dict(available_vacant_midpoint_bodies=len(C))),flush=True)
    Z=np.array(list(C),np.int32)
    tally=Counter()
    for i in range(0,len(Z),512):
        V=Z[i:i+512]
        body_dots=(V.astype(np.float32)@P.astype(np.float32).T).astype(np.int32)
        ok=body_dots.max(axis=1)<=128
        tally["pass_paired_body"]+=int(ok.sum())
        fz=frozen_constraints(V[ok])
        for b,d,bok in zip(V[ok],body_dots[ok],fz):
            k=np.flatnonzero(d>64)
            for c in np.flatnonzero(bok):
                tally["line_pass_frozen"]+=1
                on=set(g[k[pcol[k]==c]])
                off=set(g[k[pcol[k]!=c]])
                if on&off:
                    tally["inconsistent_components"]+=1
                    continue
                if on&pinned:
                    tally["frozen_component_conflict"]+=1
                    continue
                rec=dict(body2=b.tolist(),target_tail=Q[c].tolist(),equator_endpoints=C[tuple(b.tolist())],
                         flipped_pairs=int(np.isin(g,list(on)).sum()),comp=[int(x) for x in sorted(on)])
                (here/"winner_26.json").write_text(json.dumps(rec)+"\n")
                print(json.dumps(rec),flush=True)
                raise SystemExit
    print(json.dumps(dict(tally)),flush=True)
