import json,hashlib,itertools
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parent.parent/'hexagon_sparse15_20261117'
assert hashlib.sha256((P/'certificate.json').read_bytes()).hexdigest()=='5bddcfc35bb30cec8257e7b1a30775f3153fc1b796a9d12cb087b34265398eb4'
c=json.loads((P/'certificate.json').read_text());H=np.array(c['phase_generator'],dtype='int64');G=np.array(c['outer_generator'],dtype='int64')
assert H.shape==(8,15) and G.shape==(24,48)
assert np.array_equal(G[:,:24],np.eye(24,dtype=int));assert np.all((G@G.T)%3==0)
K=np.zeros((1,15),dtype=int)
for h in H:K=np.concatenate([K,(K+h)%3,(K+2*h)%3])
assert len(set(map(tuple,K)))==6561 and min(sum(x!=0) for x in K if x.any())==5
assert set(map(tuple,K))==set(map(tuple,np.load(P/'phase_code.npy')))
for r in c['profiles']:
 w=np.load(P/r['file'])['words'];assert w.shape==(r['min15_words'],r['body_blocks']);assert np.isin(w,[0,1,2]).all();assert np.all(np.sum(w!=0,axis=1)==15)
 assert len(set(map(bytes,w.astype('uint8'))))==len(w)
 f=np.zeros((len(w),48),dtype='int64');f[:,sorted(set(range(48))-set(r['shortened_coordinates']))]=w
 assert np.array_equal(f,(f[:,:24]@G)%3);assert r['sparse_points']==len(w)*6561
 print(r['body_blocks'],len(w),r['sparse_points'])
roots=[(3,1),(-3,1),(0,-2)];tri=[(2,0),(-1,1),(-1,-1)]
dot=lambda x,y:x[0]*y[0]+3*x[1]*y[1]
assert all(dot(r,r)==12 for r in roots)
assert all(abs(dot(r,t))<=6 for r in roots for t in tri)
assert all(dot((x[0]+y[0],x[1]+y[1]),(x[0]+y[0],x[1]+y[1]))>=12 for x in roots for y in roots)
assert all(dot((x[0]-y[0],x[1]-y[1]),(x[0]-y[0],x[1]-y[1]))==36 for x,y in itertools.permutations(roots,2))
assert len(set(roots+[(-a,-b) for a,b in roots]))==6
print('PASS independently reconstructed finite checks; Pless minimum15 remains explicit source premise.')
