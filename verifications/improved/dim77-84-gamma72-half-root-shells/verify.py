"""Replay with optional explicit path to dim73-82-gamma72-depth6-lattice."""
from pathlib import Path
from itertools import combinations, product
from hashlib import sha256
from collections import Counter
from sympy import Matrix
import runpy, sys, json, numpy as np

here=Path(__file__).resolve().parent
base=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else here.parent/'dim73-82-gamma72-depth6-lattice'
pins={
  'verify.py':'35b943ae846ba62a0270f1b3d5ff4f5280c208de4f4c2c5dd5df24ce3bc787df',
  'theta.py':'168d5d052ce0a6e62085800a1fa956e88791519349bab7309409dec85420f073',
  'certificate.json':'e5e656d728ef6ca1a7242e2863c8c41447247a9165a7b5aa862ea15d11a1584e',
  'results.json':'f4042a060aad10f6a3766aff8347e05fd5a9541a517516d3af3f00e49d1b5ccf',
  'comparisons.json':'d04816fd2d1cd4fa1efcc880245a3fcd104d0254c74dfaebe2b1dcf7cdeb7e4d',
  'catalog_gram.json':'1791fab8c152302db2652ec66f227eb473ce4007f934431ef6bbb2ab2f2d87ab',
}
for f,h in pins.items():
  assert sha256((base/f).read_bytes()).hexdigest()==h, f
sys.path.insert(0,str(base))
P=runpy.run_path(str(base/'verify.py'))
old=json.loads((here/'old_D4.json').read_text())
new=json.loads((here/'phases.json').read_text())
assert len(old)==len(new)==4
synd,orth,M4=P['syndrome'],P['orth'],P['M4']
for phase in old+new:
  raw=phase['raw']
  assert len(raw)==72 and set(raw)<={0,1}
  assert synd(raw)==tuple(phase['syndrome'])
H=np.array([p['raw'] for p in old+new],dtype='i8')
R=np.array([[1,-1,0,0],[0,1,-1,0],[0,0,1,-1],[0,0,1,1]],dtype='i8')
assert np.array_equal(R@R.T,P['d4'])
assert abs(Matrix(R.tolist()).det())==2 and np.all(R.sum(1)%2==0)
def coords(x):
  a,b,c,d=x
  ans=(a,a+b,(a+b+c-d)//2,(a+b+c+d)//2)
  assert tuple(np.array(ans)@R)==x
  return ans
def n10(a,b,t):
  S=(a,b,a^b)
  assert all(S)
  # Six ordered ways to place the norm-four and norm-six blocks.
  # Norm-six classes are exactly the anisotropic ones, with two representatives.
  return 4*sum(not orth(q,x^t) for p in S for q in S if p!=q for x in M4[p])
def n6(a,b,t):
  # Norm six in the Leech lattice is a single-block contribution; 4+4 is larger.
  labels=[]
  if not b: labels.append(a)
  if not a: labels.append(b)
  if a==b: labels.append(a)
  return 2*sum(not orth(i,t) for i in labels)
def lower_gate(raw,depth):
  if not depth: return
  a,b,t=synd(raw)
  S=(a,b,a^b)
  if not any(S):
    assert t and depth<=4 # A=0: no norm-four, B is singular, so no norm-six.
  elif not all(S):
    assert depth<=4
    if depth>=3: assert t not in M4[a or b]
    if depth>=4: assert orth(a or b,t)
  else:
    assert depth<=6
    if depth>=5:
      assert all(t!=x^y for p,q in combinations(S,2) for x in M4[p] for y in M4[q])
    if depth==6: P['prove'](raw,6)

# When v<=8, z lies in the radius sqrt(32/3) D4 shell, so |z_i|<=3.
# Then |x_i|<=sqrt(8)+3/2<5. The same Cartesian box includes the k=0 base.
BOX=[(x,coords(x)) for x in product(range(-4,5),repeat=4) if sum(x)%2==0]
base76=6218175600
by0=Counter()
for x,c in BOX:
  v=sum(u*u for u in x)
  if not 0<v<=8: continue
  raw=(np.array(c)@H[:4])%2
  P['verify_short'](raw,v)
  by0[v]+=P['contact'](raw,8-v)
base76+=sum(by0.values())
assert base76==P['report'][3]['count']==6221459976

counts={}
stock=Counter()
for n in range(1,5):
  patterns=Counter()
  for z,k in BOX:
    if not any(z) or any(k[n:]): continue
    h3=3*sum(u*u for u in z)
    if h3>32: continue
    for x,c in BOX:
      v4=h3+sum((2*u+w)**2 for u,w in zip(x,z))
      if v4>32: continue
      assert v4%4==0
      p=tuple(t%2 for t in c+k[:n])
      patterns[v4//4,p]+=1
  by=Counter()
  for (v,p),multiplicity in sorted(patterns.items()):
    target=8-v
    raw=(np.array(p)@H[:4+n])%2
    if n==4: stock[v,synd(raw)]+=multiplicity
    lower_gate(raw,target)
    if target==5: c=n10(*synd(raw))
    elif target==3: c=n6(*synd(raw))
    elif target==1: c=0
    else: c=P['contact'](raw,target)
    by[v]+=multiplicity*c
  total=base76+sum(by.values())
  counts[76+n]={'total':total,'hst':dict(by)}
  print(f"PASS {76+n}D kissing code {total}, gain vs pinned old result {total-P['report'][3+n]['count']}, extra contacts {dict(by)}",flush=True)
assert json.loads(json.dumps(counts))==json.loads((here/'counts77_80.json').read_text())
# Add the zero/new-tail layers of the 80D block into the complete total tail histogram.
for x,c in BOX:
  v=sum(u*u for u in x)
  if v<=8: stock[v,synd((np.array(c)@H[:4])%2)]+=1
def combine(A,B):
  U=Counter()
  for (v,(a,b,t)),k in A.items():
    for (w,(x,y,u)),j in B.items():
      if v+w<=8: U[v+w,(a^x,b^y,t^u)]+=k*j
  return U
from functools import lru_cache
def raw_of(a,b,t):
  return [w>>i&1 for w in (a^t,b,0) for i in range(24)]
@lru_cache(None)
def contact(a,b,t,w):
  raw=raw_of(a,b,t)
  assert synd(raw)==(a,b,t)
  lower_gate(raw,w)
  if w==5: return n10(a,b,t)
  if w==3: return n6(a,b,t)
  if w==1: return 0
  return P['contact'](raw,w)
def total_raw(A):
  by=Counter()
  for (v,(a,b,t)),n in A.items():
    if v:
      by[v]+=n*contact(a,b,t,8-v)
  return {'hst':dict(by),'total':6218175600+sum(by.values())}
results={80:total_raw(stock)}
assert results[80]['total']==counts[80]['total']
S=json.loads((here/'sides.json').read_text())
assert len(S)==4
for block,f in enumerate(S):
  assert len(f)==4
  for p in f:
    assert len(p['raw'])==72 and set(p['raw'])<={0,1}
    assert synd(p['raw'])==tuple(p['syndrome'])
  rows=np.array([p['raw']for p in f],dtype='i8')
  for n in range(1,5):
    A=Counter()
    for x,k in BOX:
      if any(k[n:]): continue
      w=sum(q*q for q in x)
      if w<=8: A[w,synd((np.array(k)@rows)%2)]+=1
    U=combine(stock,A)
    dim=80+4*block+n
    results[dim]=total_raw(U)
    print(f"PASS {dim}D kissing code {results[dim]['total']}, vertical contact histogram {results[dim]['hst']}",flush=True)
    if n==4: stock=U
assert json.loads(json.dumps(results))==json.loads((here/'counts80_96.json').read_text())


rows=json.loads((here/'comparisons.json').read_text())
assert [r['dimension'] for r in rows]==list(range(77,85))
a={**json.loads((here/'counts77_80.json').read_text()),**json.loads((here/'counts80_96.json').read_text())}
for r in rows:
  assert a[str(r['dimension'])]['total']==r['lower_bound']
  assert r['lower_bound']-r['preceding_repository_bound']==r['gain']>0
  assert r['lower_bound']>r['derived_comparison']
print('PASS all eight filed comparison rows')
