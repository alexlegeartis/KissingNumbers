import json
from pathlib import Path
from hashlib import sha256
import numpy as np
P=Path(__file__).resolve().parent
def rd(f,sha=None):
 b=(P/f).read_bytes()
 if sha is not None:assert sha256(b).hexdigest()==sha
 return json.loads(b)
parent=rd('pure.json','425bccd1e940053f9a34f4ea4bd24fdf1b4c3917c4756e1beaa8302d689e21f7')
weighted=rd('weighted.json','68f8687a602a51559bdbfecac75cacf90097a47e49b70bfdf711091c22527cf7')
tets=rd('tet.json','cc4151354d9b6eaf6551275a2d9424eb5b8eba11e0a4c3bf5db0bec31e9aa053')
A=np.array(parent['poles'],dtype='i8');d=rd('winner2.json')
C0=np.r_[np.delete(A,d['remove'],axis=0),np.array(d['add'],dtype='i8')]
assert C0.shape==(49,11) and np.array_equal(C0,d['poles']) and d['remove_weighted']==[]
C=np.r_[C0,-C0]
j=rd('signed_winner.json')
Y=np.r_[np.delete(C,j['removed'],axis=0),np.array(j['added'],dtype='i8')]
assert np.array_equal(Y,j['poles']) and Y.shape[1]==11 and abs(Y).max()<1000
assert len(j['removed'])==len(set(j['removed'])) and all(type(v) is int and 0<=v<len(C)for v in j['removed'])
R=j['removed_weighted']
assert len(R)==len(set(R)) and all(type(v) is int and 0<=v<32 for v in R)
N=(Y*Y).sum(1);assert np.all(N>0)
def test(P,V,M):
 q=Y@V.T
 return ((q<=0)|(4*q*q<=N[:,None]*M[None,:]))
Q=Y@Y.T;i,k=np.triu_indices(len(Y),1)
assert ((Q[i,k]<=0)|(4*Q[i,k]**2<=N[i]*N[k])).all()
assert (3*(np.abs(Y)@np.r_[np.full(8,50),np.full(3,49)])**2<=27203*N).all()
S=np.array(tets['59']['supports'],dtype='i8')
assert (5*np.abs(Y[:,S]).sum(2)**2<=6*N[:,None]).all()
Z=np.delete(np.array(weighted['weighted_tails'],dtype='i8')[:,1:],R,axis=0)
L=int(weighted['scale'])
assert 11*int(abs(Y).max())*int(abs(Z).max())<2**30
assert test(Y,Z,np.full(len(Z),L*L,dtype='i8')).all()
bound=56969157+len(Y)-94-2*len(R)
assert bound>56969157
assert rd('claims.json')['lower_bounds']=={'59':bound}
print('PASS 59D lower bound',bound,'pure signed points',len(Y),'removed weighted rows',R)
