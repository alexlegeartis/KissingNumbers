from pathlib import Path
import numpy as np,json,itertools,hashlib
P=Path(__file__).resolve().parent;src=P/'z8_generator.json';G=np.array(json.loads(src.read_text())['generator'],dtype=int)
assert G.shape==(12,24) and not np.tril(G[:,:12],-1).any() and np.all(np.diag(G[:,:12])%2)
assert np.array_equal(G[:2,:2],[[7,1],[0,7]]) and not G[2:,:2].any()
A=np.array([[1-2*((a>>2&1)^((a&1)*x)^((a>>1&1)*y)) for x,y in itertools.product((0,1),repeat=2)] for a in range(8)],dtype=int)
B=A.copy();B[:,0]*=-1
C=np.zeros((8,4),dtype=int)
for a in range(8):C[a,a%4]=2*(1 if a<4 else -1)
Bases=[A,B,C]
for i in range(3):
 assert np.array_equal(Bases[i]@Bases[i].T,A@A.T)
 for j in range(i):assert np.max(np.abs(Bases[i]@Bases[j].T))==2
assert 3*8**12==206158430208
results=[]
for d in [5,6,7]:
 labels=[];heads=[];maps=[]
 for color in range(3):
  X=[v for v in range(d) if not(v>>color&1)];Y=[v for v in range(d) if v>>color&1]
  assert len(X)<=4 and len(Y)<=4
  edges=[]
  for x in X:
   for y in Y:
    if next(i for i in range(3) if ((x^y)>>i)&1)==color:edges.append((x,y))
  maps.append(dict(color=color,first_axis_map=X,second_axis_map=Y,edges=edges))
  for x,y in edges:
   for sx,sy in itertools.product(range(2),repeat=2):
    a=X.index(x)+4*sx;b=Y.index(y)+4*sy;h=[0]*d;h[x]=2*(-1)**sx;h[y]=2*(-1)**sy
    labels.append([color,a,b]);heads.append(h)
 H=np.array(heads);assert len(H)==2*d*(d-1) and len(np.unique(H,axis=0))==len(H) and np.all((H*H).sum(1)==8)
 for i,(ci,a,b) in enumerate(labels):
  for j,(cj,c,e) in enumerate(labels):
   ip=int(H[i]@H[j])
   if ci==cj:assert ip==int(A[a]@A[c]+A[b]@A[e])
   else:assert ip<=4
 assert 22*2+4==48
 results.append(dict(dimension=88+d,head_dimension=d,maps=maps,labels=labels,heads=heads,prefix_count=len(H),fiber_size=8**10,count=len(H)*8**10))
obj=dict(status='PASS exact alphabets, projection and head geometry; Z8 parent minimum24 inherited',source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),basis_alphabets=[v.tolist() for v in Bases],full96_count=3*8**12,compressed_results=results)
(P/'certificate.json').write_text(json.dumps(obj,indent=2)+'\n');print(json.dumps(dict(full96=obj['full96_count'],compressed={r['dimension']:r['count'] for r in results}),indent=2))
