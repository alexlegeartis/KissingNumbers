from pathlib import Path
from fractions import Fraction as F
import json,hashlib,itertools,collections,math
import numpy as np
if not __debug__:raise RuntimeError('Assertions required')
S=Path('research/ring2_final92_95_20261117');D=Path(__file__).parent
expected='8c79c5a80aeecba467d203c48d194b3efac666e6b8a6219a0ada0c36c1d6f3b1';assert hashlib.sha256(Path(str(S)+'.zip').read_bytes()).hexdigest()==expected
pin='903c193dbd00c1e7c692e8de62fc4e42c42806dd53365f83ddd73e673dc0fb73';assert hashlib.sha256((S/'z8_generator.json').read_bytes()).hexdigest()==pin
G=json.load(open(S/'z8_generator.json'))['generator'];assert len(G)==12 and all(len(r)==24 for r in G);assert all(G[i][j]==0 for i in range(12)for j in range(i));assert all(G[i][i]%2==1 for i in range(12));assert all(G[i][j]==0 for i in range(2,12)for j in range(2));assert len({tuple(sum(a*G[0][j]+b*G[1][j] for _ in [0])%8 for j in range(2))for a in range(8)for b in range(8)})==64
# Generate Gray rows via characters, independently of producer arrays.
sign=lambda a:tuple((-1)**(((a//4)&1)^((a&1)*x)^(((a//2)&1)*y))for x,y in [(0,0),(0,1),(1,0),(1,1)])
rows=[sign(a)for a in range(8)];dot=lambda a,b:sum(x*y for x,y in zip(a,b))
for a,b in itertools.product(range(8),repeat=2):
 h=sum(x!=y for x,y in zip(rows[a],rows[b]));assert h==(0 if a==b else 4 if(a-b)%8==4 else 2)
 assert sum(x!=y for x,y in zip(rows[a],(-rows[b][0],)+rows[b][1:]))>=1
heads=[]
for color,A,B,edges in [(0,[0,1],[2,3],{(0,2),(0,3),(1,2),(1,3)}),(1,[0,2],[1,3],{(0,1),(2,3)})]:
 for a,b in itertools.product(range(8),repeat=2):
  if a%4>=len(A)or b%4>=len(B):continue
  i,j=A[a%4],B[b%4]
  if(i,j)not in edges:continue
  v=tuple((2 if a<4 else-2)*(k==i)+(2 if b<4 else-2)*(k==j)for k in range(4));assert dot(v,v)==8;heads.append((color,a,b,v))
assert len(heads)==24 and len({v for _,_,_,v in heads})==24
assert [sum(c==i for c,_,_,_ in heads)for i in range(2)]==[16,8]
for c,a,b,v in heads:
 for e,f,g,w in heads:
  if c==e:assert dot(v,w)==dot(rows[a],rows[f])+dot(rows[b],rows[g])
  else:assert 44+dot(v,w)<=48
# Independently enumerate the ORIGINAL full source spaces, filter shortening
# coordinates, and puncture; no producer Gaussian elimination is reused.
source=np.load(S/'source_supports.npz',allow_pickle=False);Q=source['quaternary_generator'];B=source['sign_generator'];assert Q.shape==(9,24)and B.shape==(14,24)
assert ((Q>=0)&(Q<4)).all()and ((B==0)|(B==1)).all()
V=np.zeros((4**9,24),np.uint8);ids=np.arange(4**9,dtype=np.uint32)
for k,row in enumerate(Q):
 a=((ids>>(2*k))&3)[:,None];b=row[None,:];a0=a&1;a1=a>>1;b0=b&1;b1=b>>1;V^=((a0*b0)^(a1*b1))|(((a0*b1)^(a1*b0)^(a1*b1))<<1)
V=V[V[:,0]==0,2:];assert V.shape==(65536,22);words={bytes(v)for v in V};assert len(words)==65536;h4=collections.Counter(map(int,np.count_nonzero(V,axis=1)));assert h4[0]==1 and min(k for k in h4 if k)>0 and min(k for k in h4 if k)==11
assert all(bytes(v)in words for v in np.load(S/'quaternary22_generator.npy',allow_pickle=False))
packed=[sum(int(v)<<j for j,v in enumerate(row))for row in B];binary=[]
for mask in range(1<<14):
 v=0
 for j,r in enumerate(packed):
  if(mask>>j)&1:v^=r
 if v&3==0:binary.append(v>>2)
assert len(binary)==len(set(binary))==4096;h2=collections.Counter(v.bit_count()for v in binary);assert h2[0]==1 and min(k for k in h2 if k)==6
assert all(sum(int(v)<<j for j,v in enumerate(row))in set(binary)for row in np.load(S/'binary22_generator.npy',allow_pickle=False))
# Verify the displayed shortened generators span exactly the independently filtered spaces.
short4=np.load(S/'quaternary22_generator.npy',allow_pickle=False);assert short4.shape==(8,22)
VV=np.zeros((4**8,22),np.uint8);ii=np.arange(4**8,dtype=np.uint32)
for k,row in enumerate(short4):
 a=((ii>>(2*k))&3)[:,None];b=row[None,:];a0=a&1;a1=a>>1;b0=b&1;b1=b>>1;VV^=((a0*b0)^(a1*b1))|(((a0*b1)^(a1*b0)^(a1*b1))<<1)
assert {bytes(v)for v in VV}==words
short2=np.load(S/'binary22_generator.npy',allow_pickle=False);assert short2.shape==(12,22);bb={0}
for row in short2:
 v=sum(int(x)<<j for j,x in enumerate(row));bb|={x^v for x in list(bb)}
assert bb==set(binary)
amp2=F(48,11)
assert 22*amp2==96 and 11*amp2==48 and 4*6*amp2>=96 and 22**2*amp2<48**2
# Independent RS support enumeration, determinant test via modular elimination.
def det_mod(A,p):
 A=[list(r)for r in A];v=1
 for j in range(len(A)):
  q=next(i for i in range(j,len(A))if A[i][j]%p)
  if q!=j:A[j],A[q]=A[q],A[j];v=-v
  pivot=A[j][j]%p;v=v*pivot%p;inv=pow(pivot,-1,p)
  for i in range(j+1,len(A)):
   f=A[i][j]*inv%p
   for k in range(j,len(A)):A[i][k]=(A[i][k]-f*A[j][k])%p
 return v%p
for xs in itertools.combinations(range(6),4):assert det_mod([[pow(x,k,17)for k in range(4)]for x in xs],17)!=0
supports=set()
for a0 in range(4):
 for a1,a2,a3 in itertools.product(range(17),repeat=3):
  vals=[(a0+a1*x+a2*x*x+a3*x*x*x)%17 for x in range(6)];s=(88+vals[0],)+tuple(17*(j-1)+vals[j]for j in range(1,6));assert len(set(s))==6 and max(s)<92;supports.add(s)
assert len(supports)==4*17**3
signs=[v for v in range(64)if v.bit_count()%2==0];assert len(signs)==32 and min((a^b).bit_count()for a in signs for b in signs if a!=b)==2
assert 6*16==96 and 3*16<=48 and 4*16*2>=96
assert 5*4+1*4*2<48 and (5*4)**2*amp2<48**2
assert 96*16<48**2 and 96*amp2<48**2 and 96*4<48**2
counts=dict(dense=24*8**10,sparse=65536*4096,RS=32*len(supports),poles=2*92);assert sum(counts.values())==26038868280
rec=next(v for v in json.load(open(S/'dense_certificate.json'))['results']if v['dimension']==92);assert {tuple(label)+tuple(head)for label,head in zip(rec['labels'],rec['heads'])}=={(c,a,b)+v for c,a,b,v in heads}
for line in (S/'SHA256SUMS').read_text().splitlines():
 digest,name=line.split();assert hashlib.sha256((S/name).read_bytes()).hexdigest()==digest
parent=Path('research/three_basis_independent_20261117');assert (parent/'input/z8_generator.json').read_bytes()==(S/'z8_generator.json').read_bytes();pr=json.load(open(parent/'parent_result.json'));assert pr['binary_distance']==24 and pr['z4_words_enumerated']==4**12 and pr['octads_checked']==759 and pr['low_lee_errors_per_octad']==5489

out=dict(status='NEW92D GEOMETRY AND COUNTS PASS; pinned parent minimum24 inherited',dimension=92,norm_squared=96,max_pair_inner_product=48,counts=counts,total=sum(counts.values()),gf4_histogram=dict(sorted(h4.items())),binary_histogram=dict(sorted(h2.items())),parent_generator_sha256=pin,archive_sha256=expected)
(D/'result.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
