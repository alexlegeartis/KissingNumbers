from itertools import product,combinations
from collections import Counter
from pathlib import Path
from hashlib import sha1,sha256
P=Path(__file__).resolve().parent
import json,runpy
def pin(p,h):
 b=(P.parent/p).read_bytes()
 assert sha1(f'blob {len(b)}\0'.encode()+b).hexdigest()==h,p
 return b
A='dim91-93-triality-matching-transfer/'
B=pin(A+'selection.json','0e2f71d7a6f2acaddc5a7a100e977958cb38dd8b')
I=json.loads(B)
assert len(I['append_flat_modes'])==696
R=pin(A+'supports.py','5159d6358b4da75a7bfd0156de308aae3b5942a1')
H=runpy.run_path(str(P.parent/(A+'supports.py')))['hanani'](8)
assert len(H)==385
H=tuple(frozenset(Q)for Q in H)
assert Counter(tuple(sorted(t))for h in H for t in combinations(h,3))==Counter(combinations(range(22),3))
Old=set((s,-1)for s in H)
Old|={(H[i//2],i%2) for i in I['append_flat_modes']}
F0='dim91-93-flat-parity-extensions/flat_supports.json'
p=json.loads(pin(F0,'b1112ec5baa8c6b3619f7b4215d8769c16c1c94c'))
idx=lambda n:frozenset(j for j in range(24)if n>>j&1)
F=[(idx(a),idx(u))for k in ('six','seven')for a,b in p[k] for u in (b,a^b)]
R=[F[i]for i in I['keep_old']]
assert len(R)==16 and all(len(h)==7 for h,t in R)
rw=pin('dim89-mixed-matching-rows/witness.txt','bae0e9fa8bc696ef5c23e0cc6e81fdd1be856adb')
Q=[(idx(int(v)),int(t))for v,t in map(str.split,rw.decode().splitlines())]
M=[frozenset(j-2 for j in Q[i][0])for i in I['append_matching']]
assert len(M)==916 and all(min(m)>=0 for m in M)
S0='dim91-93-complementary-parity-steiner-append/'
w=pin(S0+'witness.txt','2e45590698345e98e54e67061daf085cad887899')
W=[]
for i,a,k in (map(int,l.split())for l in w.splitlines()):
 h=idx(a)
 assert h==H[i]
 assert 0<k<15
 t=frozenset(j for b,j in enumerate(sorted(h))if k>>b&1)
 assert t and h-t
 W += [(h,t),(h,h-t)]
assert len(W)==len(set(W))==40
assert len(set(h for h,t in W))==20
for (a,p),(b,q) in combinations(W,2):
 assert sum(2 if (i in p)!=(i in q) else 4 for i in a&b)<=8
for a,p in W:
 for b,q in Old:
  cap=sum(2 if q==-1 or ((i in p)!=(q==1)) else 4 for i in a&b)
  assert cap<=8
 for b,q in R:
  cap=sum(2 if (i in p)!=(i in q) else 4 for i in a&b)
  assert cap*cap<=16*len(b)
 for b in M:
  assert 8*len(a&b)**2<=16*len(b)
all_modes=[]
for a,b,c in (map(int,l.split()) for l in (P/'witness.txt').read_text().splitlines()):
 assert all(0<=x<1<<22 for x in (a,b,c))
 assert not (a&b or a&c or b&c) and (a|b|c).bit_count()==4
 h=idx(a|b|c)
 assert h not in H
 all_modes.append({j:2 if a>>j&1 else int(bool(c>>j&1)) for j in h})
assert len(all_modes)==14 and len({tuple(sorted(x.items()))for x in all_modes})==14
capa=lambda u,v:sum(4 if u[i]==v[i] else 2 for i in u.keys()&v.keys())
for u,v in combinations(all_modes,2):assert capa(u,v)<=8
for u in all_modes:
 for b,q in Old:
  assert capa(u,{i:q if q>=0 else 2 for i in b})<=8
 for b,q in W:
  assert capa(u,{i:int(i in q) for i in b})<=8
 for b,q in R:
  assert capa(u,{i:int(i in q) for i in b})**2<=16*len(b)
 for b in M:
  assert 8*len(u.keys()&b)**2<=16*len(b)
# Norm-16 vectors against inherited dense, Golay, and ordinary roots.
assert 4*(4*4)**2<=16*88 and 4*(4*4)**2<=16*80 and 4*2**2<=16*2
# Independent exact integer row Gram matrices for flat parity frames
# and old axial, D4 doubling; choose first coordinate positive.
dot=lambda x,y:sum(i*j for i,j in zip(x,y))
T=[[v for v in product((-1,1),repeat=4)if sum(x<0 for x in v)%2==i]for i in range(2)]
Ax=[tuple(2*s if k==j else 0 for k in range(4))for j in range(4)for s in (-1,1)]
D4=[tuple(e if k==i else f if k==j else 0 for k in range(4))for i,j in combinations(range(4),2)for e,f in product((-1,1),repeat=2)]
assert all(abs(dot(x,y))<=2 for x in T[0]+T[1] for y in Ax+D4)
assert all(abs(dot(x,y))<=2 for x in Ax for y in D4)
assert all(abs(dot(x,y))==2 for x in T[0] for y in T[1])
for f in T:
 U=sorted(x for x in f if x[0]==1)
 assert len(U)==4 and all(dot(x,y)==0 for x,y in combinations(U,2))
# Construct the 1024 distinct literal mixed sign code words, and ensure
# same-frame row dot caps by direct signed symbol Gram, rather than distance.
Y=[(a,b,c,a^b^c)for a,b,c in product(range(4),repeat=3)]
Z=[tuple((v+4*(s>>i&1))for i,v in enumerate(w))for w in Y for s in range(16)]
assert len(Z)==len(set(Z))==1024
G=lambda x,y: 4*(1 if (x^y)<4 else -1)if (x&3)==(y&3)else 0
assert max(sum(G(x,y)for x,y in zip(a,b))for a,b in combinations(Z,2))==8
gain=len(Z)*len(all_modes)
assert gain==14336
prev=json.loads(pin(S0+'claims.json','4055f8326f454ca528215c0d47a1e2ef9c43200e'))
out=dict(modes=len(all_modes),points_per_mode=len(Z),
    claims=[dict(dimension=r['dimension'],previous=r['lower_bound'],gain=gain,lower_bound=r['lower_bound']+gain)for r in prev['claims']])
assert out==json.loads((P/'claims.json').read_text())
print(json.dumps(out,indent=2))
