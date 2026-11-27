import hashlib, itertools, json, math, pathlib
from fractions import Fraction
import numpy as np
from sympy import Matrix
here=pathlib.Path(__file__).resolve().parent
f=(here/"d29_certificate.npz").read_bytes()
assert hashlib.sha256(f).hexdigest()=="b27ec26a1d186ed6e9b6885921abf3b2e6bbde132206fde2689dfdaa7e39f824"
q=np.load(here/"d29_certificate.npz",allow_pickle=False)
A=q["added_integer"].astype(np.int32); B=A[:,:24]; T=A[:,24:];tail_sq=np.einsum("ij,ij->i",T,T)
pure=set(map(tuple,T[~B.any(axis=1)].tolist()))
assert all(tuple((s*8 if i==j else t*8 if i==k else 0) for i in range(5)) in pure
           for j,k in itertools.combinations(range(5),2) for s,t in itertools.product((-1,1),repeat=2))
U, ti=np.unique(T,axis=0,return_inverse=True)
M=np.flatnonzero(tail_sq==32)
_, inv, cnt=np.unique(B[M],axis=0,return_inverse=True,return_counts=True)
I=M[cnt[inv]==1]
def det(w):
    if len(w)==1:return w[0][0]
    if len(w)==2:return w[0][0]*w[1][1]-w[0][1]*w[1][0]
    if len(w)>3:
        return sum((-1)**j*a*det([x[:j]+x[j+1:] for x in w[1:]]) for j,a in enumerate(w[0]) if a)
    a,b,c=w
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
def nums(w,h):
    return tuple(det([row[:j]+(hi,)+row[j+1:] for row,hi in zip(w,h)]) for j in range(len(w)))
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def get(i,drop=None):
    d=B@B[i];d[i]=-100000
    if drop is not None:d[drop]=-100000
    if d.max()>64:return {"i":int(i),"mandatory":int(np.sum(d>64))}
    m=np.full(len(U),-100000,dtype=np.int32)
    np.maximum.at(m,ti,d)
    zero=U[m==64]
    basis=Matrix(zero.tolist() if len(zero) else [[0]*5]).nullspace()
    if not basis:return {"i":int(i),"n":0,"maxsq":"0"}
    basis=np.array([[int(x*math.lcm(*(int(u.q) for u in v))) for x in v] for v in basis],dtype=np.int64).T
    C=U@basis
    bounds={}
    for row,h in zip(C,64-m):
        if not row.any() or h >= 100000:continue
        row=tuple(map(int,row));sgn=next(x for x in row if x)
        row=tuple((x if sgn>0 else -x) for x in row)
        g=math.gcd(*row); h=Fraction(int(h),g); row=tuple(x//g for x in row)
        if row not in bounds or h<bounds[row]:bounds[row]=h
    L=math.lcm(*(h.denominator for h in bounds.values()))
    bounds=[(x,int(h*L)) for x,h in bounds.items()]
    n=basis.shape[1]; best=Fraction(0); bestvec=None
    for idx in itertools.combinations(range(len(bounds)),n):
        w,b=zip(*(bounds[i] for i in idx)); dd=det(w)
        if not dd:continue
        for ss in itertools.product((1,-1),repeat=n-1):
            bh=(b[0],)+tuple(s*x for s,x in zip(ss,b[1:]))
            p=nums(w,bh);den=abs(dd)*L
            if any(abs(dot(x,p))>h*abs(dd) for x,h in bounds):continue
            y=tuple(dot(map(int,row),p) for row in basis)
            nsq=Fraction(dot(y,y),den*den)
            if nsq>best:best,bestvec=nsq,y
    o={"i":int(i),"n":n,"planes":len(bounds),"maxsq":str(best)}
    if best>=32:
        assert bestvec is not None
        w=np.array(bestvec,dtype=np.int64); ww=int(w@w); vals=T@w
        assert all(32*int(v)**2<=ww*(64-int(k))**2 for v,k in zip(vals,d))
        o.update(exact_tail_direction=bestvec,raw_tail_squared_norm=32)
    return o
if __name__=="__main__":
    for i in I:print(json.dumps(get(i)),flush=True)
