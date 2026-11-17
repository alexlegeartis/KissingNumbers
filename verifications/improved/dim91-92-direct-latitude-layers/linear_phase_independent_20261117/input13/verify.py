from pathlib import Path
import numpy as np,json,itertools
P=Path(__file__).resolve().parent
G=np.load(P/'generator.npy').astype(int);K=np.load(P/'source_phase7812.npy');H=np.load(P/'phase13.npy')
assert np.array_equal(G[:,:24],np.eye(24,dtype=int))
assert np.array_equal(H,np.delete(K[np.all(K[:,[0,3]]==[0,2],axis=1)],[0,3],axis=1))
assert H.shape==(888,13) and len(np.unique(H,axis=0))==888
hist=np.zeros(14,dtype=np.int64)
for i in range(len(H)-1):hist+=np.bincount(np.count_nonzero(H[i+1:]!=H[i],axis=1),minlength=14)
assert not hist[:5].any() and hist[5]>0
for x,y in itertools.product((1,2),repeat=2):
 B=np.load(P/f'outer_{x}{y}0.npy');assert B.shape==(6930,45) and np.all(B<3)
 assert np.all(np.count_nonzero(B,axis=1)==13) and len(np.unique(B,axis=0))==6930
 F=np.column_stack((np.tile([x,y,0],(len(B),1)),B)).astype(int)
 assert np.array_equal((F[:,:24]@G)%3,F)
f=[0,1,-1]
for p in itertools.product(range(3),repeat=2):
 for q in itertools.product(range(3),repeat=2):
  assert sum((f[a]-f[b])**2 for a,b in zip(p,q))>=sum(a!=b for a,b in zip(p,q))
assert 4*6930*888==24615360
assert 22274253420+24615360==22298868780
assert 4*270<2025
print(json.dumps(dict(status='PASS',phase_count=888,phase_histogram=hist.tolist(),added=24615360,total92=22298868780),indent=2))
