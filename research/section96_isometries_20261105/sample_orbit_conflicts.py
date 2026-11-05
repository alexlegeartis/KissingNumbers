import numpy as np,subprocess,re,sys
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
mapping=np.fromfile(sys.argv[1],np.uint8).reshape(24,8);m=int(sys.argv[2]);mask=int(sys.argv[3]);trials=int(sys.argv[4]);s=mask.bit_count()
g=np.array([7,1,4,4,2,7,3,3,4,7,2,1]);G=np.zeros((12,24),int)
for i in range(12):G[i,i:i+12]=g;G[i,23]=3
gray=np.array([sum(((((z>>2)&1)^((z&1)*x)^(((z>>1)&1)*y))<<(2*x+y)) for x in (0,1) for y in (0,1)) for z in range(8)])
wt=np.array([int(i).bit_count() for i in range(16)])
assert mapping.shape==(24,8) and all(len(set(map(int,row)))==8 and all(int(v)<16 for v in row) for row in mapping)
assert all(wt[mapping[j,a]^mapping[j,b]]==wt[gray[a]^gray[b]] for j in range(24) for a in range(8) for b in range(8))
assert all(wt[mapping[j,a]]%2==0 for j in range(m) for a in range(8))
assert mask&((1<<m)-1)==0 and mask<1<<12
L={};R={};edges=set();witness=[]
for k in range(trials):
 p=subprocess.run(['./oracle_index',sys.argv[1],'20','0','0',str(m),str(s),str(mask),str(2026110600+k)],capture_output=True,text=True)
 if p.returncode:raise RuntimeError(p.stderr)
 q=re.search(r'WITNESS A ([0-7]{12}) B ([0-7]{12}) distance (\d+)',p.stdout)
 if q:
  a=np.array(list(q[1]),int);b=np.array(list(q[2]),int);x=a@G%8;y=b@G%8
  d=int(wt[mapping[np.arange(24),x]^gray[y]].sum());assert d==int(q[3]) and d<24
  assert all(a[i]%2==0 for i in range(12) if mask>>i&1)
  assert np.all(wt[mapping[np.arange(m),x[:m]]]==2) and np.all(wt[gray[y[:m]]]==2)
  u=tuple(a%4);v=tuple(b%4);i=L.setdefault(u,len(L));j=R.setdefault(v,len(R));edges.add((i,j));witness.append([q[1],q[2],d])
 if k%100==99:print('samples',k+1,'edges',len(edges),'candidate orbits',len(L),'source orbits',len(R),flush=True)
if edges:
 i,j=zip(*edges);A=coo_matrix((np.ones(len(edges),np.uint8),(i,j)),shape=(len(L),len(R))).tocsr();match=maximum_bipartite_matching(A,perm_type='column');print('MATCHING',int((match>=0).sum()),'exact sampled conflict edges',len(edges),flush=True)
np.savez_compressed('section92_orbit_sample.npz',mapping=mapping,witness=np.array(witness,dtype='U24'),m=m,mask=mask)
