"""Independent exact 38D audit; Python standard library only."""
import hashlib,json,itertools as I
from collections import Counter
from pathlib import Path

if not __debug__:
    raise RuntimeError("Run without Python -O")
P = Path(__file__).resolve().parent
M = P.parent / "dim37-38-qr-merge"
Q = P.parent / "dim38-qr-free-offsets"

J=json.load((P / 'witness.json').open());H=[int(v,16)for v in json.load((M / 'dense_parent39.json').open())['generators_hex']]
assert len(H)==15 and all(0<h<1<<39 for h in H)
def span(G):
 A=[0]
 for v in G:A.extend([v^a for a in A])
 return A
V=span(H);assert len(set(V))==32768
K=[v>>2 for v in V if v%4==0]; assert len(K)==8192 and min(v.bit_count()for v in K if v)==12
# Recreate the source syndrome from reduced leading-pivot combinations;
# find representatives by expanding the independent coordinate columns.
R={}
for a in V:
 if a:
  p=a.bit_length()-1
  if p not in R or a<R[p]:R[p]=a
assert len(R)==15
F=[i for i in range(39)if i not in R]
def syndrome(x):
 for p,v in sorted(R.items(),reverse=True):
  if x>>p&1:x^=v
 return sum((x>>i&1)<<j for j,i in enumerate(F))
W={}
for j in range(37):
 a=syndrome(1<<(2+j)); b=1<<j
 while a:
  p=a.bit_length()-1
  if p not in W:W[p]=a,b;break
  x,y=W[p];a^=x;b^=y
assert len(W)==24
def representative(a):
 b=0
 for p,(x,y) in sorted(W.items(),reverse=True):
  if a>>p&1:a^=x;b^=y
 assert not a
 return b
V=[]
for h,line in zip((3,0,-3),(Q / 'offset_witness.txt').open(),strict=True):
 k,*S=map(int,line.split());assert k==len(S)==len(set(S))
 for s in S:
  b=representative(s); assert syndrome(b<<2)==s
  V.append((h,b))
assert len(V)==31 and Counter(h for h,b in V)=={3:11,0:10,-3:10}
def safe(a,n,m):return a<=0 or 4*a*a<=n*m
for i,(a,c) in enumerate(V):
 for j,(b,d) in enumerate(V[:i+1]):
  dist=min((c^d^k).bit_count()for k in K if i!=j or k)
  assert dist>0 and safe(148-8*dist+a*b,148+a*a,148+b*b)
raw=(M / 'supports_3088.txt').read_bytes()
assert hashlib.sha256(raw).hexdigest()=="ccb530f90e914df84cbc8d19a71d30c79146d23a2e8e86e7a813d4d92a5a6a8b"
O=[tuple(i for i,c in enumerate(r)if c=='1')for r in raw.decode().splitlines()]
assert len(O)==3088
ix=set()
for r in J['replacements']:
 i=r['i']; assert i not in ix and tuple(r['remove'])==O[i]
 ix.add(i);O[i]=tuple(r['add'])
assert len(ix)==8
A=[set(x)for x in O]
assert len(set(map(frozenset,A)))==3088 and all(len(s)==8 and min(s)>=0 and max(s)<38 for s in A)
assert max(len(x&y)for x,y in I.combinations(A,2))==4
T=[]
for c,*l in J['stars']:
 assert len(l)==4 and len(set([c,*l]))==5 and 0<=min(c,*l)<=max(c,*l)<38
 for e in (-1,1):
  for s in I.product((-1,1),repeat=4):
   if s.count(-1)%2:continue
   T.append({c:4*e,**{j:2*v for j,v in zip(l,s)}})
assert len(T)==len({tuple(sorted(z.items()))for z in T})==272
R=[]
removed=[]
def dot(v,w):return sum(a*w.get(i,0)for i,a in v.items())
for i,j in I.combinations(range(38),2):
 for a,b in I.product((-4,4),repeat=2):
  r={i:a,j:b}
  if any(dot(t,r)>16 for t in T):removed.append((i,j,a,b))
  else:R.append(r)
assert len(R)==2584 and len(removed)==228
for a,b in I.combinations(T,2):assert dot(a,b)<=16
# Use the absolute support bound against all signed octad realizations.
for s in T:
 for o in A: assert 2*sum(abs(s.get(i,0))for i in o)<=16
# Roots/octads both norm 32; root-root distinct dot <=16 (shared one
# coordinate contributes 16; two shared sites must have at least one opposite).
# Octad parity gives <=16 on a shared support; different supports overlap <=4.
# Root-octad dot is at most 16 because at most two sites contribute 4*2.
for h in (-3,0,3):
 n=148+h*h
 for o in A:assert safe(sum(2*(abs(h)if i==0 else 2)for i in o),n,32)
 for q in R+T:assert safe(sum(abs(v)*(abs(h)if i==0 else 2)for i,v in q.items()),n,32)
total=len(V)*len(K)+128*len(O)+len(T)+len(R)
assert total==652072
print(json.dumps(dict(status='PASS',dimension=38,total=total,dense=len(V)*len(K),octads=128*len(O),new_stars=len(T),retained_roots=len(R),removed_roots=len(removed),removed_root_supports=len({(i,j) for i,j,a,b in removed})),sort_keys=True))
