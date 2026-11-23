from pathlib import Path
import json,runpy,sys,hashlib
import numpy as np
from refill import Model,sizes
assert __debug__
P=Path(__file__).resolve().parent
BASE=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else P.parent
DIMS=tuple(map(int,sys.argv[2:])) or tuple(map(int,json.loads((P/'claims.json').read_text())['lower_bounds']))
def pin(p,h):
 assert hashlib.sha256(p.read_bytes()).hexdigest()==h,(str(p),h)
src=BASE/'dim81-85-gamma72-apn-complete-dm'
pin(src/'verify.py','f9888cb194f2b16624bbabbccdedb802800b13c266a06a12710aef66fc190e33')
pin(src/'witness.json','5f89656cd4916c2534f9ed4c372ae625cf95d6fe82bfca50dd91af1b911ee334')
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
def members(Z):
 Y=set(map(tuple,Z)); assert len(Y)==len(Z) and all(tuple(-r)in Y for r in Z)
 return Z[Z[np.arange(len(Z)),np.argmax(Z!=0,axis=1)]>0]
def pos(g,d,k):
 return ((g>=0)&(d>=0)&((g>0)|(d>0)))|((g>0)&(d<0)&(g*g>k*d*d))|((g<0)&(d>0)&(k*d*d>g*g))
def tail(n):
 if n==82:
  J=np.load(P/'tail82.npz',allow_pickle=False); assert np.array_equal(J['c'],[1]*5+[3]*5)
  assert (J['B']==0).all()and J['A'].shape==(510,10) and np.max(abs(J['A']))<=144
  Z=members(J['A'].astype('i8'));g=(Z*J['c'])@Z.T;d=np.zeros_like(g);k=1;norm=144
 elif n==83:
  z=np.loadtxt(P/'tail83.txt',dtype='i8');assert z.shape==(604,22) and np.max(abs(z))<=6
  Z=members(z).reshape(-1,11,2);a,b=Z.transpose(2,0,1)
  g=a@a.T+2*(b@b.T);d=a@b.T+b@a.T;k=2;norm=36
 else:
  assert n==85
  def parse(t):return (0,int(t.split('*')[0]))if '*' in t else (int(t),0)
  Z=np.array([[parse(t)for t in r.split(',')]for r in (P/'ze99.txt').read_text().splitlines()],dtype='i8')
  assert Z.shape==(1154,13,2) and np.max(abs(Z))<=4
  Z=members(Z.reshape(len(Z),-1)).reshape(-1,13,2);a,b=Z.transpose(2,0,1)
  g=a@a.T+3*(b@b.T);d=a@b.T+b@a.T;k=3;norm=16
 assert np.array_equal(np.diag(g),np.full(len(g),norm)) and not np.diag(d).any()
 X=g-np.eye(len(g),dtype='i8')*norm
 assert not (pos(X-norm//2,d,k)|pos(-X-norm//2,-d,k)).any()
 E=pos(g,d,k);F=pos(-g,-d,k)
 return len(g),np.block([[E,F],[F,E]])
for n in DIMS:
 nn,EG=tail(n)
 Z=np.load(P/f'witness{n}.npz',allow_pickle=False)
 rows=[tuple(map(int,r))for r in Z['row']]
 assert len(rows)==2*nn and {i for i,c,t in rows}==set(range(2*nn))
 assert all(0<=c<3 and 0<=t<128 for i,c,t in rows)
 s,m,v=M.N
 assert Z['triple'].shape==(2*nn,(s+7)//8) and Z['full66'].shape==(2*nn,(v+7)//8) and Z['full48'].shape==(2*nn,m)
 assert Z['triple'].dtype==Z['full66'].dtype==np.dtype('u1') and Z['full48'].dtype==np.dtype('u4')
 A=np.unpackbits(Z['triple'],axis=1,bitorder='little',count=s).astype(bool)
 B=Z['full48'];C=np.unpackbits(Z['full66'],axis=1,bitorder='little',count=v).astype(bool)
 assert not (B>>24).any()
 occ=list(zip(A,B,C));tot=np.zeros(3,dtype='i8')
 for j,(i,c,t) in enumerate(rows):
  u=[k for k,(v,p,q) in enumerate(rows[:j])if abs(i-v)!=nn and p==c and (q==t or EG[i,v])]
  X=M.add(i,c,t,[rows[k]for k in u],[occ[k]for k in u],EG)
  assert all(not (x&~y).any()for x,y in zip(occ[j],X)),(n,j)
  tot+=sizes(occ[j])
 count=6218175600+2*nn+int(tot@[8,4,4])
 J=json.loads((P/'claims.json').read_text())
 assert count==J['lower_bounds'][str(n)] and count>J['previous'][str(n)]
 print('PASS exact signed-tail mask replay',n,count,'profile counts',tot.tolist(),'gain',count-J['previous'][str(n)],flush=True)
