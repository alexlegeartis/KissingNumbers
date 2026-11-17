from pathlib import Path
import json,itertools,hashlib
import numpy as np
P=Path(__file__).resolve().parent;src=P/'z8_generator.json';G=np.array(json.loads(src.read_text())['generator'],dtype=np.int64)
assert G.shape==(12,24) and not np.tril(G[:,:12],-1).any() and np.all(np.diag(G[:,:12])%2)
assert np.array_equal(G[:2,:2],[[7,1],[0,7]]) and not G[2:,:2].any()
assert len({tuple((np.array([a,b])@G[:2,:2])%8) for a,b in itertools.product(range(8),repeat=2)})==64
Gray=np.array([[(a>>2&1)^((a&1)*x)^((a>>1&1)*y) for x,y in itertools.product((0,1),repeat=2)] for a in range(8)],dtype=np.int16)
S=1-2*Gray;F=S.copy();F[:,0]*=-1
assert np.all(Gray.sum(1)%2==0) and np.all(((1-F)//2).sum(1)%2==1)
assert np.min(np.count_nonzero(S[:,None,:]!=F[None,:,:],axis=2))==1
assert np.array_equal(S@S.T,F@F.T)
models={7:[([0,1,2],[3,4,5,6],None),([0,3,4],[1,2,5,6],[(0,1),(0,2),(3,5),(3,6),(4,5),(4,6)])],4:[([0,1],[2,3],None),([0,2],[1,3],[(0,1),(2,3)])],5:[([0,1],[2,3,4],None),([0,2],[1,3,4],[(0,1),(2,3),(2,4)])],6:[([0,1,2],[3,4,5],None),([0,3],[1,2,4,5],[(0,1),(0,2),(3,4),(3,5)])]}
results=[]
for d,modelspec in models.items():
 labels=[];heads=[]
 for coset,(A,B,edges) in enumerate(modelspec):
  assert set(A).isdisjoint(B) and len(A)<=4 and len(B)<=4
  for a,b in itertools.product(range(8),repeat=2):
   ax,ay=a%4,b%4
   if ax>=len(A) or ay>=len(B):continue
   x,y=A[ax],B[ay]
   if edges is not None and (x,y) not in edges:continue
   h=[0]*d;h[x]=2*(1 if a<4 else -1);h[y]=2*(1 if b<4 else -1)
   labels.append([coset,a,b]);heads.append(h)
 H=np.array(heads);assert np.all((H*H).sum(1)==8) and len(np.unique(H,axis=0))==len(H)
 for i,(ci,a,b) in enumerate(labels):
  for j,(cj,c,e) in enumerate(labels):
   ip=int(H[i]@H[j])
   if ci==cj:assert ip==int(S[a]@S[c]+S[b]@S[e])
   else:assert 16-2*ip>=8
 counts=[sum(p[0]==i for p in labels) for i in (0,1)];expected={7:([48,24],77309411328),4:([16,8],25769803776),5:([24,12],38654705664),6:([36,16],55834574848)}
 assert counts==expected[d][0] and len(labels)*8**10==expected[d][1]
 results.append(dict(dimension=88+d,head_dimension=d,maps=modelspec,labels=labels,heads=heads,prefix_counts=counts,fiber_size=8**10,count=len(labels)*8**10))
obj=dict(status='PASS exact finite maps, Gram and projection checks; parent minimum24 inherited',source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),norm_squared=96,results=results)
(P/'dense_certificate.json').write_text(json.dumps(obj,indent=2)+'\n');print(json.dumps({r['dimension']:dict(count=r['count'],prefix_counts=r['prefix_counts']) for r in results},indent=2))
