#!/usr/bin/env python3
"""Retain every previous point and share some Γ72 phases with new 85D tail rays."""
from pathlib import Path
from hashlib import sha256
import sys,os,runpy,json,tempfile
import numpy as np
assert __debug__
here=Path(__file__).resolve().parent
parent=here.parent/'dim82-85-bw16-tail-sections'
pins={'verify.py': '28dc2b85d228045e6974e5ede2d940ff09c2beebf23be8ac99a18e3bdbaa7975', 'phases.json': '050157ff3f8ee6282e4092dd208bbd8ca1f83ef4b4a6efdda9bdc16735b1f21f', 'claims.json': 'a7ecc17a84cea1ca130e50b75f5273a5ac5ac95eb788569e0ac05b394eba5f1a'}
assert all(sha256((parent/n).read_bytes()).hexdigest()==s for n,s in pins.items())
with tempfile.TemporaryDirectory() as tmp:
 os.environ['EXPORT_TRIPLE']=tmp
 R=runpy.run_path(str(parent/'verify.py'))
 O=np.load(Path(tmp)/'triple384.npy').astype(bool)
 assert O.shape==(384,13824)
os.environ.pop('EXPORT_TRIPLE')
assert R['out']['85']==6292777568 and int(O.sum())==4028344
rows=R['rows'];T=R['ts'];H=R['H'];rb=R['rb'];V=R['tails'](13)
assert len(V)==453
G=V@V.T
assert np.all(np.diag(G)==8) and np.abs(G-8*np.eye(len(V),dtype='i8')).max()<=4
W=R['phases']['85']
S=[[i,c,t]for i,c,t in W];occ=list(O)
assert len(S)==len(occ)==384 and len({(c,t)for i,c,t in S})==len(S)
assert all(tuple(V[i])==tuple(v)and c+1==d and H[t]==u for (i,c,t),(v,d,u) in zip(S,rows))
K={};o={}
for t,u in enumerate(H):
 A,B=T[rb[u]]
 y=(A[:,2].astype('u8')<<24)|(B[:,0].astype('u8')<<12)|B[:,1].astype('u8')
 assert len(np.unique(y))==len(y)
 K[t]=y;q=np.argsort(y);o[t]=y[q],q
extra=json.loads((here/'extra.json').read_text())
add=0
for i,c,t,count in extra:
 assert 0<=i<len(V) and i not in {j for j,d,u in S}
 assert (c,t) in {(d,u)for j,d,u in S[:384]}
 A=np.ones(13824,bool)
 for k,(j,d,u) in enumerate(S):
  if d!=c:continue
  if t==u:
   A&=~occ[k] # no two tails share the same exact horizontal vector
  elif G[i,j]!=0:
   idx=np.flatnonzero(A)
   if not len(idx):break
   delta=rb[H[t]^H[u]];assert delta
   Q=K[t][idx][None,:]^np.array([0,delta,delta<<12],dtype='u8')[:,None]
   keys,p=o[u]
   J=np.minimum(np.searchsorted(keys,Q),len(keys)-1)
   A[idx]&=~((keys[J]==Q)&occ[k][p[J]]).any(axis=0)
 assert int(A.sum())==count>0,(i,count,int(A.sum()))
 add+=count
 S.append((i,c,t));occ.append(A)
total=6292777568+16*add+2*(len(V)-384)
claims={'lower_bounds':{'85':total}}
assert claims==json.loads((here/'claims.json').read_text())
print('PASS shared phases',total,'new triples',add,'new tilted rays',len(extra),'new pure contacts',2*(len(V)-384))
