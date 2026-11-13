import os
os.environ['OPENBLAS_NUM_THREADS']='2'
import numpy as np, json, hashlib, time
from pathlib import Path
p=Path(__file__).resolve().parent
data=p.parent/'data'
def h(name):return hashlib.sha256((data/name).read_bytes()).hexdigest()
x=np.load(data/'coordinates30_248180.npz')['points']; z=np.load(data/'delta30_248180.npz'); a=z['added'];r=z['removed'];s=np.load(data/'source30_248084.npz')['points']
assert np.issubdtype(x.dtype,np.integer) and x.shape==(248180,30)
assert s.shape==(248084,30)
for q in [x,s,a,r]:
 assert np.issubdtype(q.dtype,np.integer) and np.all(q>=-8) and np.all(q<=8)
assert not np.any(x%2) and not np.any(s%2)
x=x//2;s=s//2
assert np.all(np.sum(s.astype(np.int64)**2,axis=1)==32)
rows=lambda q:set(map(tuple,q.tolist()))
X,S,A,R=map(rows,[x,s,a,r]);assert len(X)==len(x) and len(S)==len(s) and len(A)==len(a) and len(R)==len(r)
assert R<=S and not A&S and X==(S-R)|A
assert len(a)==384 and len(r)==288
assert np.all(np.sum(x.astype(np.int64)**2,axis=1)==32)
# Exact binary64 products: coordinates integral; each product integer,
# Cauchy-Schwarz bounds sum of absolute products by32. All partial sums
# and products are therefore representable, no tolerance or rounding premise.
f=x.astype(np.float64);n=len(x);pairs=0;mx=-32;t=time.time();B=2048
for i in range(0,n,B):
 q=f[i:i+B];d=q@q.T;np.fill_diagonal(d,-32);assert d.max(initial=-32)<=16
 pairs+=len(q)*(len(q)-1)//2;mx=max(mx,int(d.max(initial=-32)))
 for j in range(i+B,n,B):
  d=q@f[j:j+B].T;assert d.max(initial=-32)<=16;pairs+=len(q)*len(f[j:j+B]);mx=max(mx,int(d.max(initial=-32)))
 if i%(B*10)==0:print(i,pairs,time.time()-t,flush=True)
assert pairs==n*(n-1)//2
# independent int64 changed-pair replay against retained source and each other
ret=np.array(sorted(S-R),dtype=np.int64);ai=a.astype(np.int64)
for i in range(0,len(ret),4096):assert (ai@ret[i:i+4096].T).max(initial=-32)<=16
d=ai@ai.T;np.fill_diagonal(d,-32);assert d.max(initial=-32)<=16
report={'pass':True,'dimension':30,'count':n,'norm_squared':32,'max_dot':mx,'unordered_pairs':pairs,'added':len(a),'removed':len(r),'source_count':len(s),'arithmetic':'binary64 exactly representable integer products and partial sums bounded by32; no tolerance. Independent int64 changed-pair replay.','hashes':{k:h(k) for k in ['coordinates30_248180.npz','delta30_248180.npz','source30_248084.npz']},'seconds':time.time()-t}
(p/'audit_result.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
