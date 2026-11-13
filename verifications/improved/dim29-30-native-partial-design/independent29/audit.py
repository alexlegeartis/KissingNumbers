import numpy as np,json,hashlib,time
from pathlib import Path
p=Path(__file__).resolve().parent;data=p.parent/'data';start=time.time();f=data/'coordinates29_231814.npz';z=np.load(f,allow_pickle=False);X=z['points'];assert X.dtype==np.int16 and X.shape==(231814,29);assert int(z['squared_norm'])==128 and int(z['inner_product_limit'])==64 and int(z['coordinate_denominator'])==2;assert np.all(np.sum(X.astype(np.int64)**2,axis=1)==128);assert len(np.unique(X,axis=0))==len(X)
s=np.load(data/'source29_231776.npz')['points'];d=np.load(data/'delta29_231814.npz');added=d['added'].astype(np.int16)*2;removed=d['removed'].astype(np.int16)*2;assert added.shape==(492,29) and removed.shape==(454,29);assert s.shape==(231776,29)
def keys(x):return {row.tobytes() for row in np.asarray(x,dtype='<i2')}
S=keys(s);R=keys(removed);A=keys(added);F=keys(X);assert len(S)==len(s) and len(R)==len(removed) and len(A)==len(added);assert R<=S;assert not(A&S);assert (S-R)|A==F;assert len(F)==231776-454+492==231814
# Binary64 products and ALL possible partial sums exact: integer inputs, d*maxabs²<2^53.
maxabs=int(np.max(np.abs(X)));bound=29*maxabs**2;assert bound<2**53;Y=X.astype(np.float64);N=len(X);pairs=0;maximum=-1e9
for k in range(0,N,128):
 e=min(k+128,N);G=Y[k:e]@Y[k:].T
 for j in range(e-k):G[j,:j+1]=-np.inf
 m=float(np.max(G));assert m<=64,(k,m);maximum=max(maximum,m);pairs+=sum(N-i-1 for i in range(k,e))
 if k%8192==0:print(k,pairs,time.time()-start,flush=True)
assert pairs==N*(N-1)//2
r=dict(pass_=True,count=N,dimension=29,squared_norm=128,inner_product_limit=64,max_inner_product=int(maximum),distinct=True,all_unordered_pairs=pairs,arithmetic='binary64 exact integer dot products; no antipodal grouping',max_partial_sum_bound=bound,integer_exact_limit=2**53,source_count=len(s),added=len(added),removed=len(removed),gain=len(added)-len(removed),coordinate_sha256=hashlib.sha256(f.read_bytes()).hexdigest(),delta_sha256=hashlib.sha256((data/'delta29_231814.npz').read_bytes()).hexdigest(),source_sha256=hashlib.sha256((data/'source29_231776.npz').read_bytes()).hexdigest(),seconds=time.time()-start)
json.dump(r,open(p/'receipt.json','w'),indent=2);print(json.dumps(r),flush=True)
