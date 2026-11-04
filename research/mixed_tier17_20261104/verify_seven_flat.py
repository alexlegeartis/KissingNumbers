import numpy as np,itertools
R='research/mixed_tier17_20261104/';V=np.fromfile(R+'mixed17_orbits.bin',np.int8).reshape(-1,16).astype(np.int16);H=np.fromfile(R+'mixed17_group.bin',np.uint16);F=np.fromfile(R+'mixed17_flat_reps.bin',np.uint16);Sg=np.array([[(-1 if int(h)>>i&1 else 1) for i in range(16)] for h in H],np.int16);words=np.array([[(-1 if (int(f)^int(h))>>i&1 else 1) for i in range(16)] for f in F for h in H],np.int16);same=(V@words.T).reshape(1920,64,32).max(axis=2)>8
flat=np.zeros((64,64),bool)
for i in range(64):
 for j in range(i):flat[i,j]=flat[j,i]=min((int(F[i])^int(F[j])^int(h)).bit_count() for h in H)<6
nb=[j for j in range(1,64) if not flat[0,j]];count=0
for x in itertools.combinations(nb,6):
 if any(flat[a,b] for a,b in itertools.combinations(x,2)):continue
 S=(0,)+x;I=np.flatnonzero(np.all(~same[:,S],axis=1));assert len(I)==16 and all(i>=896 for i in I)
 for a,b in itertools.combinations(I,2):assert int(np.max(Sg@(V[a]*V[b])))>8
 count+=1
assert count==56;print('PASS: all 56 seven-flat-coset sets containing zero leave 16 mutually conflicting tier orbits; at most one can be chosen.')
