"""Independent exact 60D and 61D shell/literal inequalities."""
from pathlib import Path
from itertools import combinations,product
from functools import reduce
from operator import xor
import hashlib,sys,subprocess,tempfile,struct
assert __debug__
parent=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent.parent/'dim62-67-rm2-selfsparse'
P=parent/'witness.txt'
assert hashlib.sha256(P.read_bytes()).hexdigest()=='5173adcf040a2c2978522f781aef88cab6cbfad7c7bcde6f1c99584554f32c8f'
F=parent/'gray.cpp'
assert hashlib.sha256(F.read_bytes()).hexdigest()=='6d08f0a74d9cca4f0697775ad38959a3f669acd3261c11377a787da366c45f1b'
with tempfile.TemporaryDirectory() as tmp:
 exe=Path(tmp)/'gray'; raw=Path(tmp)/'supports.bin'
 subprocess.run(['g++','-std=c++17','-O3','-march=native',F,'-o',exe],check=True)
 subprocess.run([exe,P,raw],check=True)
 S={x for x, in struct.iter_unpack('=Q',raw.read_bytes())}
assert len(S)==37996 and all(x.bit_count()==16 for x in S)
C=list(map(int,P.read_text().split()))
M=[sum(1<<x for x in range(64) if x&i==i) for i in range(64)]
L=[v for i,v in enumerate(M) if i.bit_count()<=2]
H=[v for i,v in enumerate(M) if i.bit_count()==3]
assert len(L)==22 and len(H)==20 and len(C)==80
W=[reduce(xor,(v for i,v in enumerate(H) if p>>i&1),0) for p in C]
assert all(not (x&y).bit_count()%2 for x in W+L for y in W+L)
assert all(not (x&y).bit_count()%2 for x in S for y in W+L)
marked=(0,20,41,61);mask=sum(1<<i for i in marked)
plus,minus=0,20
J=[x for x in S if not x&mask]
A=[x for x in S if x&mask==1<<plus]
B=[x for x in S if x&mask==1<<minus]
assert (len(J),len(A),len(B))==(12005,3956,3956)
Z={0}
for v in L:
 s=sum(((v>>i)&1)<<j for j,i in enumerate(marked));Z|={k^s for k in Z}
assert len(Z)==16
# Odd/uniform 16- and 15-coordinate extended Hamming sign banks.
for n in (15,16):
 D=[[]for _ in range(32)]
 for a in range(1<<n):
  s=reduce(xor,(16+j for j in range(n) if a>>j&1),0)
  D[s].append(a)
 assert [len(x)for x in D]==[1<<(n-5)]*32
 assert min(x.bit_count()for x in D[0]if x)==4
 assert all(x.bit_count()%2==j//16 for j,a in enumerate(D)for x in a)
# Unit-sphere compatibility is tested with squared norms and exact max dot.
def safe(dot,p,q):
 assert dot<=0 or 4*dot*dot<=p*q,(dot,p,q)
# 60D dense; linearly cut rank4 in each coset, any two labels differ <=1.
for x,y in product((0,1),repeat=2):safe(28+2*(x^y).bit_count(),60,60)
safe(16-2,60,16) # even versus odd parity forces at least one disagreement
safe(4*4,60,64)
safe(8,60,64)
# 61D dense vectors scaled by 2; labels 0..3, vertical heights 1,4,-1,-4.
hs=[1,4,-1,-4]
for x,h in enumerate(hs):
 for y,k in enumerate(hs):
  safe(4*(28+2*(x^y).bit_count())+h*k,240+h*h,240+k*k)
 safe(2*(16-2),240+h*h,16)
 for t,m in [(1,plus),(-1,minus)]:
  rank=marked.index(m)
  risky=1 if t>0 else 3
  forced=1-((risky>>rank)&1)
  for y,k in enumerate(hs):
   bit=(y>>rank)&1
   disag=2 if bit!=forced else 0
   safe(2*(15-disag)+k*t,240+k*k,16)
 safe(max(4*4*2,4*(3*2+abs(h))),240+h*h,64)
 safe(8*max(2,abs(h)),240+h*h,64)
# Dense support types are distinct from the other shells; q4 affine plane
# stocks have ≤2 common sites, otherwise identical support has ≥one sign flip.
for d,p in [(8,16),(16,64),(32,64)]:
 safe(d,p,64 if p!=16 else 16)
safe(8,16,16);safe(16,16,64);safe(32,64,64)
# Spectral supports and the exact carry planes, independently counted.
P4={sum(1<<x for x in (a,b,c,a^b^c)) for a,b,c in combinations(range(64),3)}
assert len(P4)==10416 and all(x.bit_count()==4 for x in P4)
p60=sum(not x&mask for x in P4)
r61=mask^1
p61=sum(not x&r61 for x in P4)
assert p60==7995 and p61==8555
n60=2*80*(1<<18)+len(J)*2048+16*p60+120
n61=4*80*(1<<18)+len(J)*2048+(len(A)+len(B))*1024+16*p61+122
assert (n60,n61)==(66657320,116711210)
print('PASS exact 60 and 61:',n60,n61,'sparse',len(J),len(A),len(B),'planes',p60,p61)
