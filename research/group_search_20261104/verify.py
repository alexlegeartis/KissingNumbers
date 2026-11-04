"""Independent exact certificates. No floating-point acceptance criterion."""
from fractions import Fraction as F
from itertools import product, combinations
from pathlib import Path
import json,hashlib
import numpy as np
import sympy as sp
P=Path(__file__).resolve().parent
# Q(sqrt(-3)); tuples a,b mean a+b sqrt(-3).
zero=(F(0),F(0));one=(F(1),F(0));omega=(F(-1,2),F(1,2))
def add(x,y):return (x[0]+y[0],x[1]+y[1])
def mul(x,y):return (x[0]*y[0]-3*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def scale(x,c):return (c*x[0],c*x[1])
def powers(x,n):
 r=[one]
 for _ in range(n-1):r.append(mul(r[-1],x))
 return r

def build10():
    pts=list(product(range(3),repeat=2));ix={x:i for i,x in enumerate(pts)};w=powers(omega,3)
    def transform(v,g):
        if g in (0,2):
            a=g//2; out=[]
            for y in pts:
                val=zero
                for b in range(3):
                    x=list(y);x[a]=b
                    val=add(val,mul(w[y[a]*b%3],v[ix[tuple(x)]]))
                out.append(mul((F(0),F(-1,3)),val))
            return tuple(out)
        a=(g-1)//2
        return tuple(mul(w[(x[a]**2 if g<4 else x[0]*x[1])%3],v[i]) for i,x in enumerate(pts))
    def orb(v):
        vs=[tuple(v)];seen={tuple(v)}
        for v in vs:
            for g in range(5):
                u=transform(v,g)
                if u not in seen:seen.add(u);vs.append(u)
            assert len(vs)<=1000
        return vs
    a=[zero]*9;a[0]=one
    b=[zero]*9
    for k in (3,6):b[k]=(F(1,2),F(0))
    for k in (4,8):b[k]=(F(1,4),F(1,4))
    A=orb(a);B=orb(b);assert (len(A),len(B))==(80,270)
    X=[tuple(mul(z,x) for x in v) for z in w for v in A]+B
    assert len(set(X))==510
    # Re(v*conj(w)) = a*c+3*b*d. Real embedding is [a,sqrt(3)b].
    # Common denominator 12 yields integer form diag(I,3I).
    den=1
    from math import lcm
    for v in X:
        for a,b in v:den=lcm(den,a.denominator,b.denominator)
    C=np.array([[int(a*den) for a,b in v]+[int(b*den) for a,b in v] for v in X],dtype=np.int64)
    return C,np.array([1]*9+[3]*9),den*den,{'orbit_sizes':[len(A),len(B)],'construction':'3 phases of 80 + 270; exact Fourier/phase action in Q(sqrt(-3))'}

def build14():
    # E7 roots scaled by 2. Complements of Fano lines {i,i+1,i+3}.
    roots=[2*np.eye(7,dtype=int)[i] for i in range(7)]
    for i in range(7):
        S=sorted(set(range(7))-{i,(i+1)%7,(i+3)%7})
        for signs in product((-1,1),repeat=3):
            v=np.zeros(7,dtype=int);v[S]=[1,*signs];roots.append(v)
    R=np.array(roots);assert len(R)==63
    adj=[set(np.flatnonzero(R@v==0)) for v in R]
    cliques=[]
    def rec(part,cand):
        if len(part)==7:cliques.append(tuple(part));return
        while len(part)+len(cand)>=7:
            a=min(cand);cand.remove(a);rec(part+[a],cand&adj[a])
    rec([],set(range(63)))
    by=[[] for _ in R]
    for c in cliques:
        for i in c:by[i].append(c)
    def cover(left,chosen):
        if not left:return chosen
        i=min(left,key=lambda i:sum(set(c)<=left for c in by[i]))
        for c in by[i]:
            if set(c)<=left:
                result=cover(left-set(c),chosen+[c])
                if result is not None:return result
    frames=cover(set(range(63)),[]);assert len(frames)==9
    X=[]
    # All ordered distinct axes in each frame, with independent signs.
    for frame in frames:
        for i in frame:
            for j in frame:
                if i!=j:
                    for s,t in product((-1,1),repeat=2):X.append(np.r_[s*R[i],t*R[j]])
    assert len(X)==1512
    for v in R:
        for s,t in product((-1,1),repeat=2):X.append(np.r_[s*v,t*v])
    for side in (0,7):
        for i,j in combinations(range(7),2):
            for s,t in product((-1,1),repeat=2):
                v=np.zeros(14,dtype=int);v[side+i]=2*s;v[side+j]=2*t;X.append(v)
    return np.array(X),np.ones(14,dtype=int),8,{'frame_count':9,'available_orthogonal_frames':len(cliques),'components':[1512,252,168]}

def check(C,weights,norm,expected_n,expected_rank):
    assert C.dtype.kind in 'iu' and weights.dtype.kind in 'iu'
    assert len(C)==expected_n and len(np.unique(C,axis=0))==expected_n
    G=(C*weights)@C.T
    assert np.all(np.diag(G)==norm)
    upper=G[np.triu_indices(len(C),1)]
    assert int(upper.max())*2<=norm
    # Exact rank of the integer rational coefficient matrix is also the rank
    # after multiplying columns by nonzero sqrt(weights).
    rank=int(sp.Matrix(C.tolist()).rank())
    assert rank==expected_rank
    return {'cardinality':len(C),'rank_exact':rank,'norm_squared_integer_form':int(norm),
            'max_inner_product_integer_form':int(upper.max()),'pairs_checked':len(upper),
            'violations':int(np.count_nonzero(2*upper>norm)),
            'sha256_integer_coordinates':hashlib.sha256(C.astype('<i8').tobytes()).hexdigest()}

def main():
    report={}
    for dim,builder,N in [(10,build10,510),(14,build14,1932)]:
        C,w,n,meta=builder();r=check(C,w,n,N,dim);r.update(meta)
        np.savez_compressed(P/f'certified_{dim}.npz',coordinates=C,weights=w,norm_squared=n)
        # Exercise rejection of duplicates and a norm-preserving sign mutation
        # that creates a violating pair.
        bad=C.copy();bad[0]=bad[1]
        try:check(bad,w,n,N,dim)
        except AssertionError:r['duplicate_rejected']=True
        else:raise AssertionError('duplicate accepted')
        report[str(dim)]=r;print(dim,r,flush=True)
    (P/'exact_verification.json').write_text(json.dumps(report,indent=2))
if __name__=='__main__':main()
