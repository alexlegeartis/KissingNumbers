"""Replay the exact 71D block refill; requires NumPy and SymPy."""
from pathlib import Path
from hashlib import sha256, sha1
from itertools import combinations
import numpy as np
import json,runpy,sys
HERE=Path(__file__).resolve().parent
P=Path(sys.argv[1]) if len(sys.argv)>1 else HERE.parent/'dim73-82-gamma72-depth6-lattice'
old=Path(sys.argv[2]) if len(sys.argv)>2 else HERE.parent/'dim68-91-lattice-code-continuation'/'mixed_shell71_20261117'/'verify71.py'
pins={
  'verify.py':'35b943ae846ba62a0270f1b3d5ff4f5280c208de4f4c2c5dd5df24ce3bc787df',
  'theta.py':'168d5d052ce0a6e62085800a1fa956e88791519349bab7309409dec85420f073',
  'certificate.json':'e5e656d728ef6ca1a7242e2863c8c41447247a9165a7b5aa862ea15d11a1584e',
  'results.json':'f4042a060aad10f6a3766aff8347e05fd5a9541a517516d3af3f00e49d1b5ccf',
  'comparisons.json':'d04816fd2d1cd4fa1efcc880245a3fcd104d0254c74dfaebe2b1dcf7cdeb7e4d',
  'catalog_gram.json':'1791fab8c152302db2652ec66f227eb473ce4007f934431ef6bbb2ab2f2d87ab',
}
for f,h in pins.items(): assert sha256((P/f).read_bytes()).hexdigest()==h, f
d=old.read_bytes(); assert sha1(f'blob {len(d)}\0'.encode()+d).hexdigest()=='263a5b243690d6480e8ea0bf7aca19d8045673da'
previous=runpy.run_path(str(old))
assert previous['total']==2623578936
assert [r[:2] for r in previous['slices']]==[(8,3),(8,4),(10,4),(10,5)]
sys.path.insert(0,str(P))
D=runpy.run_path(str(P/'verify.py'))
data=json.loads((HERE/'pattern.json').read_text())
l=np.array(data['lambda'],dtype='i8'); W=D['W']; F=D['F']; A=D['A']; B=D['B']; M4=D['M4']
assert l.shape==(24,) and sum(np.all(W==l,axis=1))==1 and l@F@l==4
dots=(W@F)@l
assert np.count_nonzero(dots==0)==data['orthogonal_minimum_leech_vectors']==93150
removed=int(np.count_nonzero(dots==2)); assert removed==data['removed']==4600
Z=W[dots==0]+l
target=((Z@A%2)@D['p24']).tolist()
phase=((Z@B%2)@D['p24']).tolist()
# Test xor two-sums by witness lookup, not a precomputed sumset.
assert len(M4)==4095 and 0 not in M4 and all(len(v)==24 for v in M4.values())
allowed=sum(bool(a) and all((b^p) not in M4[a] for p in M4[a]) for a,b in zip(target,phase))
assert allowed==data['added']==81952
# Confirm all possible first-block squared norms of a raw norm-sixteen triple.
legal={0,4,6,8,10,12,14,16}
triples=[(a,b,c) for a in legal for b in legal for c in legal if a+b+c==16
         and (0 not in (a,b,c) or sorted((a,b,c)) in ([0,0,16],[0,8,8]))]
assert set(map(lambda t:tuple(sorted(t)),triples))=={(0,0,16),(0,8,8),(4,4,8),(4,6,6)}
assert all(t[0] not in (10,12,14) for t in triples)
total=previous['total']+allowed-removed
assert total==data['lower_bound']==2623656288 and data['prior']==previous['total']
print('PASS exact 71D lower bound',total,'add',allowed,'remove',removed)
