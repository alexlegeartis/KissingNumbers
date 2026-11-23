#!/usr/bin/env python3
"""Exact 82/83D Γ72 tail transfer."""
from pathlib import Path
import hashlib,json,runpy,sys
import numpy as np
from refill import Model,sizes
from components import improve

assert __debug__
P=Path(__file__).resolve().parent
BASE=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else P.parent
def pin(p,h):
 b=p.read_bytes()
 assert hashlib.sha256(b).hexdigest()==h,(str(p),h)
 return b
src=BASE/'dim81-85-gamma72-apn-complete-dm'
pin(src/'verify.py','f9888cb194f2b16624bbabbccdedb802800b13c266a06a12710aef66fc190e33')
pin(src/'witness.json','5f89656cd4916c2534f9ed4c372ae625cf95d6fe82bfca50dd91af1b911ee334')
previous=json.loads(pin(BASE/'dim82-85-bw16-tail-sections'/'claims.json',
    'a7ecc17a84cea1ca130e50b75f5273a5ac5ac95eb788569e0ac05b394eba5f1a'))['lower_bounds']
# The published exact replay supplies literal norm 4, 6 and 8 Leech
# classes and all horizontal phase layouts. No sample stock is used.
sys.argv=[str(src/'verify.py')]
R=runpy.run_path(sys.argv[0])
H=R['H'];AA=np.array(sorted(R['RA'],key=R['RA'].get),dtype='u4')
BB=np.array(sorted(R['RB'],key=R['RB'].get),dtype='u4')
Q=[]
for t in H:
 A,B=R['T'][R['RB'][t]]
 Q.append(np.stack((AA[A[:,2]],BB[B[:,0]],BB[B[:,1]]),axis=1))
Q=np.array(Q)
U=np.array([R['layouts'][t][0] for t in H])
U[:,:,0]=AA[U[:,:,0]]
U[:,:,1]=BB[U[:,:,1]]
V=np.array([R['layouts'][t][1] for t in H])
M=Model(H,*[R[k]for k in ('A8','A6','F','V8','V6')],U,V,Q)

def reps(Z):
 V=set(map(tuple,Z))
 assert len(V)==len(Z) and all(tuple(-r)in V for r in Z)
 L=[i for i,r in enumerate(Z) if r[np.flatnonzero(r)[0]]>0]
 assert 2*len(L)==len(Z)
 return L

def pos(a,b):
 # Exact positivity in Z[sqrt(2)].
 if a==0:return b>0
 if b==0:return a>0
 if a>0:return b>0 or a*a>2*b*b
 return b>0 and 2*b*b>a*a

def tail(d):
 if d==82:
  p=P/'tail82.npz'
  pin(p,'bc4d17399173208c664700e5a97ff6e8b04ffe54b72b5ce20a54c951e9e8015c')
  T=np.load(p)
  assert np.all(T['B']==0) and T['c'].tolist()==[1]*5+[3]*5
  Z=T['A'].astype('i8')
  assert Z.shape==(510,10) and np.max(np.abs(Z))<=144
  q=Z[reps(Z)]
  G=(q*T['c'])@q.T
  assert np.all(np.diag(G)==144) and np.max(np.abs(G-144*np.eye(255,dtype='i8')))<=72
  E=G!=0
 else:
  Z=np.loadtxt(P/'tail83.txt',dtype='i8')
  assert Z.shape==(604,22) and np.max(np.abs(Z))<=6
  q,p=Z[reps(Z)].reshape(302,11,2).transpose(2,0,1)
  A=q@q.T+2*(p@p.T); B=q@p.T+p@q.T
  for i in range(302):
   assert A[i,i]==36 and B[i,i]==0
   for j in range(i):
    a,b=map(int,(A[i,j],B[i,j]))
    assert not pos(a-18,b) and not pos(-a-18,-b)
  E=(A!=0)|(B!=0)
 np.fill_diagonal(E,False)
 return E

out={}
for d in (82,83):
 E=tail(d);N=len(E)
 rows=json.loads((P/f'assign{d}.json').read_text())
 assert len(rows)==N and {i for i,c,t in rows}==set(range(N))
 assert len({(c,t)for i,c,t in rows})==N
 assert all(0<=c<3 and 0<=t<128 for i,c,t in rows)
 print('PASS exact tail geometry',d,2*N,'nonorthogonal ray pairs',int(E.sum()//2),flush=True)
 pre=[];occ=[];tot=np.zeros(3,dtype='i8')
 for i,c,t in rows:
  z=M.add(i,c,t,pre,occ,E)
  pre.append((i,c,t));occ.append(z);tot+=sizes(z)
 k=improve(M,Q,E,rows,occ)
 count=6218175600+2*N+16*(int(tot[0])+k)+8*int(tot[1]+tot[2])
 assert count>previous[str(d)]
 print('PASS',d,count,'raw profile totals',tot.tolist(),'reallocated groups',k,'delta',count-previous[str(d)],flush=True)
 out[str(d)]=count
assert out==json.loads((P/'claims.json').read_text())['lower_bounds']
print(json.dumps({'status':'PASS','lower_bounds':out},sort_keys=True))
