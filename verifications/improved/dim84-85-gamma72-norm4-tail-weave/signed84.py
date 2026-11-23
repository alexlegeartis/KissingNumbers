from pathlib import Path
import json,runpy,sys,hashlib
import numpy as np
from refill import Model,sizes
assert __debug__
P=Path(__file__).resolve().parent
BASE=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else P.parent
DIMS=() if len(sys.argv)>2 and sys.argv[2]=='--geometry' else (tuple(map(int,sys.argv[2:])) or tuple(map(int,json.loads((P/'base84.json').read_text())['lower_bounds'])))
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
 elif n==84:
  Z=np.array(json.loads((P/'tail84.json').read_text()),dtype='i8')
  assert Z.shape==(420,12) and np.max(abs(Z))<=2
  assert len(np.unique(np.concatenate((Z,-Z)),axis=0))==840
  g=Z@Z.T;d=np.zeros_like(g);k=1;norm=4
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
 J=json.loads((P/'base84.json').read_text())
 assert count==J['lower_bounds'][str(n)] and count>J['previous'][str(n)]
 print('PASS exact signed-tail mask replay',n,count,'profile counts',tot.tolist(),'gain',count-J['previous'][str(n)],flush=True)


from functools import lru_cache
from refill import pack24
@lru_cache(9000)
def signs(i,j):
 D=M.R8[i]@M.V8[j].T
 return pack24(D==4),pack24(D==-4)
@lru_cache(None)
def rel(a,b):
 d=M.delta[a,b];assert d
 a0,b0=M.match(1,a,b,d)
 a1,b1=M.match(1,a,b,np.uint64(0))
 assert np.all(M.U[a,a0,2]==M.U[b,b0,2])
 assert np.all(M.U[a,a1,:2]==M.U[b,b1,:2])
 tab=np.array([signs(int(i),int(j))for i,j in zip(M.U[a,a1,2],M.U[b,b1,2])],dtype='u4').reshape(-1,2,24)
 return a0,b0,a1,b1,tab
base=json.loads((P/'base84.json').read_text())
new=json.loads((P/'apron84.json').read_text())
for n,c in new['lower_bounds'].items():
 nn,E=tail(int(n));Z=np.load(P/f'witness{n}.npz',allow_pickle=False)
 S=np.load(P/f'apron{n}.npz',allow_pickle=False)
 A=S['full48'];assert A.dtype==np.dtype('u4') and A.shape==(len(Z['row']),M.N[1],2)
 assert not (A>>24).any() and np.array_equal(S['row'],Z['row'])
 rows=np.array(Z['row'],dtype='i8')
 for j,(i,v,a) in enumerate(rows):
  m=A[j]
  for k,(h,w,b) in enumerate(rows[:j]):
   if w!=v or abs(i-h)==nn:continue
   if a==b:assert not (m&A[k]).any()
   elif E[i,h]:
    a0,b0,a1,b1,Sign=rel(a,b)
    assert not (m[a0]&A[k,b0]).any()
    if not len(a1):continue
    X0=Sign[:,0];X1=Sign[:,1]
    P0=A[k,b1,0,None];P1=A[k,b1,1,None]
    V0=pack24(((X0&P0)|(X1&P1))!=0)
    V1=pack24(((X1&P0)|(X0&P1))!=0)
    assert not (V0&m[a1,0]).any() and not (V1&m[a1,1]).any()
 before=sum(int(x).bit_count()for x in Z['full48'].ravel())
 after=sum(int(x).bit_count()for x in A.ravel())
 gain=2*after-4*before
 assert gain>0 and c==base['lower_bounds'][n]+gain and c>base['previous'][n]
 print('PASS signed 4+8 replacement',n,c,'gain over fixed signed',gain,flush=True)
