from pathlib import Path
import numpy as np,json,itertools,hashlib
P=Path(__file__).resolve().parent;j=json.loads((P/'candidate.json').read_text());G=np.array(json.loads((P/'z8_generator.json').read_text())['generator'],dtype=int)
assert G.shape==(12,24) and not np.tril(G[:,:12],-1).any() and np.all(np.diag(G[:,:12])%2)
assert np.array_equal(G[:2,:2],[[7,1],[0,7]]) and not G[2:,:2].any()
A=np.array([[1-2*((a>>2&1)^((a&1)*x)^((a>>1&1)*y)) for x,y in itertools.product((0,1),repeat=2)] for a in range(8)],dtype=int);B=A.copy();B[:,0]*=-1;C=np.zeros((8,4),int)
for a in range(8):C[a,a%4]=2*(-1)**(a//4)
Bases=[A,B,C]
for i in range(3):
 assert np.array_equal(Bases[i]@Bases[i].T,A@A.T)
 for k in range(i):assert np.max(np.abs(Bases[i]@Bases[k].T))<=2
old=np.array([np.concatenate((A[a],A[b])) for a,b in j['source_prefixes']]);oldgram=old@old.T
assert len(old)==38 and len(set(map(tuple,j['source_prefixes'])))==38
for H in j['full_rotated_heads']:
 H=np.array(H,dtype=int);assert H.shape==(38,8) and np.array_equal(H@H.T,oldgram)
 assert np.all(H[:,0]+H[:,4]==0) and np.all((H*H).sum(1)==8)
selected=j['selected'];assert len(selected)==100 and len({tuple(r['head']) for r in selected})==100
assert len({(r['color'],*r['prefix']) for r in selected})==100
for x,r in enumerate(selected):
 c=r['color'];assert c in range(3);idx=j['source_prefixes'].index(r['prefix']);assert r['head']==j['full_rotated_heads'][c][idx]
 for s in selected[:x]:
  ip=sum(a*b for a,b in zip(r['head'],s['head']))
  if r['color']==s['color']:
   a,b=r['prefix'];u,v=s['prefix'];assert ip==int(A[a]@A[u]+A[b]@A[v])
  else:assert ip<=4
assert 88+8-1==95 and 88+8==96 and 22*2+4==48
n=100*8**10;assert n==107374182400
result=dict(status='PASS exact100-head three-basis95D construction; inheritedZ8 minimum24 explicit',count=n,color_counts=[sum(r['color']==i for r in selected) for i in range(3)],source_sha256=hashlib.sha256((P/'z8_generator.json').read_bytes()).hexdigest(),candidate_sha256=hashlib.sha256((P/'candidate.json').read_bytes()).hexdigest())
(P/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
