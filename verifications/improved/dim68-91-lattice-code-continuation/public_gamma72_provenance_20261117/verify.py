"""Exact public Gamma72 source and norm8-completion verifier. Requires NumPy.
Optionally pass a downloaded neb72.html to verify the source matrix extraction.
"""
import hashlib,json,re,sys
from pathlib import Path
from math import comb
import numpy as np
P=Path(__file__).resolve().parent;s=json.load(open(P/'source.json'));G,F,A,B,T=[np.array(s[k],dtype=np.int64)for k in ['G','F','A','B','T']];I=np.eye(24,dtype=np.int64);Z=np.zeros((24,24),dtype=np.int64)
canon=('\n'.join(' '.join(map(str,r))for r in G)+'\n').encode();assert hashlib.sha256(canon).hexdigest()==s['gram72_canonical_ascii_sha256']
if len(sys.argv)>1:
 raw=Path(sys.argv[1]).read_bytes();part=raw.decode().split('<a NAME="GRAM">')[1].split('<a NAME="DET">')[0];v=list(map(int,re.findall(r'-?\d+',re.sub('<[^>]*>',' ',part))));assert v[0]==72 and np.array_equal(np.array(v[1:]).reshape(72,72),G);print('Source HTML matrix matches; raw hash:',hashlib.sha256(raw).hexdigest())
assert np.array_equal(T,np.block([[I,I,B],[Z,A,A],[Z,Z,2*I]]));D=np.kron(np.eye(3,dtype=np.int64),F);assert np.array_equal(T@D@T.T,2*G)
assert np.array_equal(A+B,I)and np.array_equal(A@B,2*I)and np.array_equal(B@A,2*I)
assert np.array_equal(A@F@A.T,2*F)and np.array_equal(B@F@B.T,2*F)and np.array_equal(A@F+F@A.T,F)
def det(M):
 M=[list(map(int,r))for r in M];n=len(M);prev=1;sign=1
 for k in range(n-1):
  if M[k][k]==0:
   j=next(j for j in range(k+1,n)if M[j][k]);M[k],M[j]=M[j],M[k];sign=-sign
  p=M[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):num=M[i][j]*p-M[i][k]*M[k][j];assert num%prev==0;M[i][j]=num//prev
  for i in range(k+1,n):M[i][k]=0
  prev=p
 return sign*M[-1][-1]
assert det(F)==1 and np.all(np.diag(F)%2==0)and np.array_equal(F,F.T)
# Theta coefficients, q exponent = squared norm /2.
E=[1]+[240*sum(d**3 for d in range(1,n+1)if n%d==0)for n in range(1,5)]
def mul(a,b):return [sum(a[i]*b[n-i]for i in range(n+1))for n in range(5)]
theta=[x-720*y for x,y in zip(mul(mul(E,E),E),[0,1,-24,252,-1472])];assert theta==[1,0,196560,16773120,398034000]
assert 1+196560//2+16773120//2+398034000//48==2**24
V=np.load(P/'leech_min4_witnesses.npz')['vectors'].astype(np.int64);assert V.shape==(196560,24)and len(np.unique(V,axis=0))==196560;assert np.all(np.sum((V@F)*V,axis=1)==4)
keys0=((V&1)*(1<<np.arange(24))).sum(1).astype(np.uint32);keys,counts=np.unique(keys0,return_counts=True);assert len(keys)==98280 and np.all(counts==2)
bkeys=((B&1)*(1<<np.arange(24))).sum(1).astype(np.uint32);groups=np.zeros(len(keys),np.uint32)
for i in range(24):groups^=((keys>>i)&1)*bkeys[i]
ug,gc=np.unique(groups,return_counts=True);assert len(ug)==4095 and np.all(gc==24)and 0 not in ug
span={0};basis=[]
for row in A:
 k=sum((int(x)&1)<<i for i,x in enumerate(row))
 if k not in span:basis.append(k);span|={x^k for x in list(span)}
assert len(span)==4096 and len(basis)==12
U=np.array(sorted(span),dtype=np.uint32);bits=(U[:,None]>>np.arange(24))&1;assert np.all((bits@B)%2==0);assert np.all((np.sum((bits@F)*bits,axis=1)//2)%2==0)
# Each nonzero B fiber is a nontrivial affine linear quadratic restriction.
reps=np.array([keys[np.flatnonzero(groups==g)[0]]for g in ug]);R=(reps[:,None]>>np.arange(24))&1;Ub=(np.array(basis)[:,None]>>np.arange(24))&1;functions=R@F@Ub.T%2;assert np.all(np.any(functions,axis=1))
# A nonzero linear functional on F2^12 has exactly2048 zeros and2048 ones.
N=np.zeros(1<<24,dtype=np.uint16)
for g in ug:
 k=keys[groups==g];np.add.at(N,(k[:,None]^k[None,:]^g).ravel(),1)
N*=4
for u in span:
 if u:N[u]+=96
assert N[0]==0 and np.all(N[keys]==0)and int(N.sum())==9828000
support=np.flatnonzero(N).astype(np.uint32);q=np.zeros(len(support),np.uint32)
def parity(x):
 x=x.copy();x^=x>>16;x^=x>>8;x^=x>>4;x^=x>>2;x^=x>>1;return x&1
for i in range(24):
 mask=np.uint32(sum((int(F[i,j])&1)<<j for j in range(i+1,24)));q^=((support>>i)&1)*(parity(support&mask)^((int(F[i,i])//2)&1))
assert not np.any(q)
assert hashlib.sha256(N.astype('<u2').tobytes()).hexdigest()==s['expected_fiber8_table_little_endian_uint16_sha256']
weights={'c_norm0':4095*48**2+2*196560,'c_norm4':48*96+2048*4,'c_norm6':48*4,'fiber8_table_sum':int(N.sum()),'c_norm16_pure':1};assert weights['c_norm0']==9828000 and weights['c_norm4']==12800 and weights['c_norm6']==192
shell8=weights['c_norm0']+196560*weights['c_norm4']+16773120*weights['c_norm6']+48*int(N.sum())+196560;assert shell8==6218175600
report={'status':'PASS','source_gram_entries_matched':5184,'gram72_canonical_ascii_sha256':s['gram72_canonical_ascii_sha256'],'Leech_theta_coefficients_through_norm8':theta,'norm4_distinct_integer_witnesses':196560,'type2_parity_classes':98280,'nonzero_B_fibers':4095,'type2_classes_per_nonzero_B_fiber':24,'quadratic_split_per_nonzero_B_fiber':[2048,2048],'norm8_completion_weights':weights,'Gamma72_norm8_total':shell8,'table_sha256':s['expected_fiber8_table_little_endian_uint16_sha256'],'premises':s['external_theorem_premises']};json.dump(report,open(P/'verification.json','w'),indent=2);print(json.dumps(report,indent=2))
