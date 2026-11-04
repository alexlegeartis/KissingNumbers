"""Independent exact audit of the 5-flat-coset clique-cover certificate."""
import itertools,json,numpy as np
R='research/mixed_tier17_20261104/';V=np.fromfile(R+'mixed17_orbits.bin',np.int8).reshape(-1,16).astype(np.int16);H=np.fromfile(R+'mixed17_group.bin',np.uint16);F=np.fromfile(R+'mixed17_flat_reps.bin',np.uint16);cert=json.load(open(R+'mixed17_five_flat_clique_certificate.json'));assert len(cert)==280 and len(set(map(int,H)))==32 and len(set(map(int,F)))==64
sign=np.array([[(-1 if int(h)>>i&1 else 1) for i in range(16)] for h in H],np.int16);words=np.array([[(-1 if (int(f)^int(h))>>i&1 else 1) for i in range(16)] for f in F for h in H],np.int16)
same=(V@words.T).reshape(1920,64,32).max(axis=2)>8
flat=np.zeros((64,64),bool)
for i in range(64):
 for j in range(i):flat[i,j]=flat[j,i]=min((int(F[i])^int(F[j])^int(h)).bit_count() for h in H)<6
assert np.array_equal(same,np.fromfile(R+'mixed17_same.bin',np.uint8).reshape(1920,64).astype(bool))
assert np.array_equal(flat,np.fromfile(R+'mixed17_flat_graph.bin',np.uint8).reshape(64,64).astype(bool))
nb=[j for j in range(1,64) if not flat[0,j]];expected={(0,)+x for x in itertools.combinations(nb,4) if not any(flat[a,b] for a,b in itertools.combinations(x,2))};assert len(expected)==280 and {tuple(x['flat']) for x in cert}==expected
for item in cert:
 S=tuple(item['flat']);allowed=set(np.flatnonzero(np.all(~same[:,S],axis=1)).tolist());cliques=item['cliques'];assert len(cliques)==16 and all(len(c)==6 for c in cliques)
 used=[int(v) for c in cliques for v in c];assert len(used)==len(set(used))==96 and set(used)==allowed
 for c in cliques:
  for a,b in itertools.combinations(c,2):
   assert int(np.max(sign@(V[a]*V[b])))>8
print('PASS: all 280 flat five-coset sets containing zero; 16 independently checked six-cliques cover the 96 eligible tier orbits in every case.')
