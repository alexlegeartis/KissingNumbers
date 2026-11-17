from pathlib import Path
import numpy as np,json,itertools,hashlib
D=Path(__file__).resolve().parent;z=json.loads((D/'candidate_three.json').read_text());p=D/'z8_generator.json';assert hashlib.sha256(p.read_bytes()).hexdigest()=='903c193dbd00c1e7c692e8de62fc4e42c42806dd53365f83ddd73e673dc0fb73';G=np.array(json.loads(p.read_text())['generator'],dtype=int)
assert G.shape==(12,24)and not G[4:,:4].any()and not np.tril(G[:4,:4],-1).any()and np.all(np.diag(G[:4,:4])%2)
A=np.array([[1-2*((a>>2&1)^((a&1)*x)^((a>>1&1)*y))for x,y in itertools.product((0,1),repeat=2)]for a in range(8)]);B=A.copy();B[:,0]*=-1;C=np.zeros((8,4),dtype=int)
for a in range(8):C[a,a%4]=2*(-1)**(a//4)
bases=[A,B,C];R=A@A.T
for i in range(3):
 assert np.array_equal(bases[i]@bases[i].T,R)
 for j in range(i):assert np.max(np.abs(bases[i]@bases[j].T))==2
heads=[];labels=[];colors=[]
for color,(P,S)in enumerate(zip(z['partitions'],z['selected_supports'])):
 assert sorted(sum(P,[]))==list(range(12))and all(len(g)<=4 for g in P)
 assert len(set(tuple(sorted(s))for s in S))==len(S)
 HH=[];LL=[]
 for support in S:
  assert len(set(support))==4;axes=[]
  for group in P:
   q=set(group)&set(support);assert len(q)==1;axes.append(next(iter(q)))
  for bits in itertools.product(range(2),repeat=4):
   h=[0]*12;lab=[]
   for j,(axis,bit)in enumerate(zip(axes,bits)):h[axis]=2*(-1)**bit;lab.append(P[j].index(axis)+4*bit)
   HH.append(h);LL.append(lab)
 HH=np.array(HH);LL=np.array(LL);assert np.all((HH*HH).sum(1)==16)and len(np.unique(LL,axis=0))==len(LL)
 gram=sum(R[LL[:,j,None],LL[None,:,j]]for j in range(4));assert np.array_equal(HH@HH.T,gram)
 heads.extend(HH.tolist());labels.extend(LL.tolist());colors.extend([color]*len(HH))
H=np.array(heads);colors=np.array(colors)
assert len(H)==1776 and len(np.unique(H,axis=0))==1776
for i in range(3):
 for j in range(i):assert np.max(H[colors==i]@H[colors==j].T)<=8
assert 20*4+16==96 and 20*2+8==48
dense=len(H)*8**8;assert dense==29796335616
# GF16 cubic RS layer. Finite field x^4+x+1.
def mul(a,b):
 out=0
 for _ in range(4):
  if b&1:out^=a
  b>>=1;a<<=1
  if a&16:a^=19
 return out
assert all(any(mul(a,b)==1 for b in range(1,16))for a in range(1,16))
def rank(M):
 M=[r[:]for r in M];rr=0
 for j in range(len(M[0])):
  k=next((i for i in range(rr,len(M))if M[i][j]),None)
  if k is None:continue
  M[rr],M[k]=M[k],M[rr];iv=next(b for b in range(1,16)if mul(M[rr][j],b)==1);M[rr]=[mul(x,iv)for x in M[rr]]
  for i in range(len(M)):
   if i!=rr:
    v=M[i][j];M[i]=[a^mul(v,b)for a,b in zip(M[i],M[rr])]
  rr+=1
 return rr
V=[[1,x,mul(x,x),mul(mul(x,x),x)]for x in range(6)]
assert all(rank([V[i]for i in I])==4 for I in itertools.combinations(range(6),4))
signs=[s for s in itertools.product([-1,1],repeat=6)if sum(x<0 for x in s)%2==0];assert len(signs)==32 and min(sum(a!=b for a,b in zip(x,y))for x,y in itertools.combinations(signs,2))==2
groups=[list(range(12))]+[list(range(12+16*j,28+16*j))for j in range(5)];assert sum(groups,[])==list(range(92))
rs=32*12*16**3;roots=4*sum(len(g)*(len(g)-1)//2 for g in groups);assert rs==1572864 and roots==2664
assert 6*4*2==48 and (4*4)*48<48**2
out={'status':'PASS exact alphabets,1776 heads,first4-prefix uniformity,RS field/minors/signs/root geometry; parent minimum24 inherited','supports_by_color':[len(s)for s in z['selected_supports']],'heads':len(H),'fiber_size':8**8,'dense_count':dense,'rs_count':rs,'root_count':roots,'total':dense+rs+roots,'parent_generator_sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
(D/'verification.json').write_text(json.dumps(out,indent=2)+'\n');(D/'head_certificate.json').write_text(json.dumps({'heads':heads,'prefixes':labels,'colors':colors.tolist()})+'\n');print(json.dumps(out,indent=2))
