"""Replay the literal append against a pinned exact half-root Gamma72 shell."""
from pathlib import Path
from functools import lru_cache
from itertools import combinations, permutations
from hashlib import sha256
import json,sys,runpy
import numpy as np
here=Path(__file__).resolve().parent
up=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else here.parent/'dim77-84-gamma72-half-root-shells'
base=Path(sys.argv[2]).resolve() if len(sys.argv)>2 else up.parent/'dim73-82-gamma72-depth6-lattice'
assert __debug__
pins={
 'verify.py':'3bcfc646f39054912043b60340c0109dd14d2c299fdd4f5e217735727566b045',
 'old_D4.json':'3eb7abf15afd4bef1b287bac44782cf6b7ac7aa9434a0d2b6054fd15d3ea2643',
 'phases.json':'41b12d3da3904373de62a05bf4466873e195bbc0bad7b1193dc1ba85accc84a8',
 'sides.json':'4730c18663c23789ee70b1b02004b78770c1cb9521312dd220f72471c85a0772',
 'counts77_80.json':'ecd3df41a8d99c2bfbe4ec0a0ea3a9d6bed775e211a62b98388279aecb805810',
 'counts80_96.json':'c334633f45d2043cd364e935dfa3170de810cd4a283c1f12873c6e75bfe3123a',
 'comparisons.json':'03f031542a6ff1343ccebcce55575677525aa99a15dbfe794be0c42f8569b585',
}
assert all(sha256((up/k).read_bytes()).hexdigest()==v for k,v in pins.items())
sys.argv=[str(up/'verify.py'),str(base)]
U=runpy.run_path(sys.argv[0])
P=U['P']
M4=P['M4']; sy=P['syndrome']; orth=P['orth']; BOX,H=U['BOX'],U['H']; H=H.astype('i8')
assert (H.shape==(8,72))
@lru_cache(None)
def depth(a,b,t):
 if not(a or b):return 4 if t else 0
 if not(a and b and a!=b):
  z=a or b
  return 2 if t in M4[z] else 4 if orth(z,t) else 3
 S=(a,b,a^b)
 if any(x^t in M4[y] for x0,y in combinations(S,2) for x in M4[x0]):return 4
 if any(not orth(y,x^t) for x0,y in permutations(S,2) for x in M4[x0]):return 5
 # Any independent pair has a two-occupied-block representation with norm <=4+8.
 return 6
def dot(x,y):return sum(a*b for a,b in zip(x[:4],y[:4]))+3*sum(a*b for a,b in zip(x[4:],y[4:]))
def F(z):
 a,b,t=z
 out=[(w>>i)&1 for w in (a^t,b,0) for i in range(24)]
 assert sy(out)==z
 return out
def xor(x,y):return tuple(a^b for a,b in zip(x,y))
def pairs_raw(n):
 for z,k in BOX:
  if any(k[n:]):continue
  sz=3*sum(a*a for a in z)
  if sz>32:continue
  for x,c in BOX:
   u=tuple(2*a+b for a,b in zip(x,z))
   p=sz+sum(a*a for a in u)
   if p>32:continue
   assert p%4==0
   coeff=np.array(c+k[:n],dtype='i8')
   phase=sy(coeff@H[:4+n]%2)
   yield p//4,phase,u+z
@lru_cache(None)
def pairs(n):return tuple(pairs_raw(n))
@lru_cache(None)
def stock_count(v,phase):
 if v==0:return 0
 w=8-v
 raw=F(phase);U['lower_gate'](raw,w)
 if w==5:return U['n10'](*phase)
 if w==3:return U['n6'](*phase)
 if w==1:return 0
 return P['contact'](raw,w)
@lru_cache(None)
def construct_stock(d):
 n=min(4,d-76);m=max(0,d-80)
 E=list(pairs(n))
 L=[]
 A=np.array([x['raw']for x in U['S'][0]],dtype='i8')
 for u,y in BOX:
  if any(y[m:]):continue
  w=sum(t*t for t in u)
  if w>8:continue
  phase=sy(np.array(y,dtype='i8')@A%2)
  for s,h,z in E:
   v=s+w
   if not 0<v<=8:continue
   p=xor(phase,h); k=stock_count(v,p)
   if k:L.append((v,p,z,k))
 return L
def check_old(v,p,x,w,r,y):
 # Old vector's vertical first coordinates are (y[:4],sqrt(3)y[4:])/2.
 # Appends have no component in the orthogonal sidecar.
 s=depth(*xor(p,r))+v+w-8
 return s>=0 and s*s*dot(x,x)>=v*dot(x,y)**2
def check_new(v,p,x,w,r,y):
 s=depth(*xor(p,r))+v+w-8
 return s>=0 and s*s*dot(x,x)*dot(y,y)>=4*v*w*dot(x,y)**2
W=json.loads((here/'witness.json').read_text())
claims=json.loads((here/'claims.json').read_text())
for dstr,rows in sorted(W.items()):
 d=int(dstr);assert 77<=d<=84
 stock=construct_stock(d)
 old=6218175600+sum(k for v,p,z,k in stock)
 published=U['counts'][d]['total']if d<80 else U['results'][d]['total']
 assert old==published
 total=old
 for i,(v,p,x,k) in enumerate(rows):
  assert len(x)==8 and all(type(z)==int for z in x) and dot(x,x)
  if d==77:assert x[4]==-x[5] and x[6]==x[7]==0
  if d==78:assert sum(x[4:7])==0 and x[7]==0
  if d==79:assert sum(x[4:])==0
  p=tuple(p)
  assert v in (4,8)
  if v==8:assert p==(0,0,0) and k==1
  if v==4:assert depth(*p)==4 and P['contact'](F(p),4)==k>0
  for w,r,y,n in stock:assert check_old(v,p,x,w,r,y)
  for w,r,y,n in rows[:i]:assert check_new(v,p,x,w,tuple(r),y)
  # Its opposite uses the identical phase because the quotient has exponent two;
  # vertical distance to it is sqrt(4*v) >= sqrt(8).
  total+=2*k
 assert claims[dstr]=={'prior':old,'gain':total-old,'lower_bound':total}
 print('PASS',d,total,'+%d'%(total-old),len(rows),'antipodal tail fibre pairs')
