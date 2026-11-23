#!/usr/bin/env python3
"""Exact replay of the 85D Zinoviev–Ericson / Gamma72 splice."""
import json, os, runpy, sys, tempfile
from hashlib import sha256
from pathlib import Path
import numpy as np
from refill import Model,sizes
from components import improve
assert __debug__
here=Path(__file__).resolve().parent
repo=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else here.parent
parent=repo/'dim82-85-bw16-tail-sections'
pins={'verify.py':'28dc2b85d228045e6974e5ede2d940ff09c2beebf23be8ac99a18e3bdbaa7975',
'phases.json':'050157ff3f8ee6282e4092dd208bbd8ca1f83ef4b4a6efdda9bdc16735b1f21f',
'claims.json':'a7ecc17a84cea1ca130e50b75f5273a5ac5ac95eb788569e0ac05b394eba5f1a'}
assert all(sha256((parent/n).read_bytes()).hexdigest()==h for n,h in pins.items())
helper=repo/'dim81-85-gamma72-half-sign/two_profiles.py'
assert sha256(helper.read_bytes()).hexdigest()=='0159dc533e090920cc3371b4744437e9371d5868ecbe3dcab103e29a88ca70d4'
with tempfile.TemporaryDirectory() as aux:
 os.environ['EXPORT_PROFILE']=aux
 sys.argv=[str(parent/'verify.py')]
 R=runpy.run_path(sys.argv[0])
 X={k:v for k,v in np.load(Path(aux)/'profile_stock.npz').items()}
 O={k:v for k,v in np.load(Path(aux)/'profile_85.npz').items()}
 del os.environ['EXPORT_PROFILE']
assert R['out']['85']==6292777568
R0=R['R'];Ainv=np.array(sorted(R0['RA'],key=R0['RA'].get),dtype='u4')
Binv=np.array(sorted(R0['RB'],key=R0['RB'].get),dtype='u4')
Q=[]
for t in X['H']:
 A,B=R0['T'][R0['RB'][int(t)]]
 Q.append(np.stack((Ainv[A[:,2]],Binv[B[:,0]],Binv[B[:,1]]),axis=1))
Q=np.array(Q)
M=Model(*[X[k]for k in ('H','A8','A6','F','V8','V6','U','V')],Q)
# Coordinates are exact pairs in Z[sqrt(3)], with squared norm 16.
Z=[]
for r in (here/'ze99.txt').read_text().splitlines():
 vals=[]
 for s in r.split(','):
  if s.endswith('*sqrt(3)'):vals.append((0,int(s.split('*')[0])))
  else:vals.append((int(s),0))
 assert len(vals)==13;Z.append(vals)
Z=np.array(Z,dtype='i8');assert Z.shape==(1154,13,2) and np.max(np.abs(Z))<=4
flat=Z.reshape(1154,26);A={tuple(r)for r in flat}
assert len(A)==1154 and all(tuple(-r)in A for r in flat)
keep=[i for i,r in enumerate(flat) if r[np.flatnonzero(r)[0]]>0]
assert len(keep)==577
a,b=Z[keep].transpose(2,0,1)
G=a@a.T+3*(b@b.T);D=a@b.T+b@a.T
def pos(g,d):
 if d==0:return g>0
 if g==0:return d>0
 if g>0 and d>0:return True
 if g<0 and d<0:return False
 return (g*g>3*d*d) if g>0 else (3*d*d>g*g)
for i in range(577):
 assert G[i,i]==16 and D[i,i]==0
 for j in range(i):
  g,d=map(int,(G[i,j],D[i,j]));assert not (pos(g-8,d) or pos(-g-8,-d))
E=(G!=0)|(D!=0);np.fill_diagonal(E,False)
rows=json.loads((here/'assign.json').read_text())
assert len({i for i,c,t in rows})==len(rows)
assert len({(c,t)for _,c,t in rows[:384]})==384
tot=np.zeros(3,dtype='i8')
P=[];occ=[]
for i,c,t in rows:
 assert 0<=i<577 and 0<=c<3 and 0<=t<128
 z=M.add(i,c,t,P,occ,E);P.append((i,c,t));occ.append(z);tot+=sizes(z)
more=improve(M,Q,E,rows,occ)
count=int(6218175600+1154+16*(tot[0]+more)+8*(tot[1]+tot[2]))
claims=json.loads((here/'claims.json').read_text())
assert count==claims['lower_bounds']['85'] and count>6295120802
print('PASS',count,'raw profile totals',tot.tolist(),'reallocated triple groups',more,'assigned rays',len(rows))
