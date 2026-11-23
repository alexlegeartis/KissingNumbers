"""Independent literal 36D Cartesian expansion, exact 64-bit integer Gram products."""
from pathlib import Path
import hashlib,json,importlib.util,sys,itertools
from collections import Counter
import numpy as np
P=Path(__file__).resolve().parent
pr=P.parent/'dim36-ct-triangle-refill'
db=P.parent/'dim36-ct-double-trine'
nc=P.parent/'dim36-ct-nonlocal'
assert __debug__
def read(p):
 return json.loads(p.read_text())
def H(p,h):
 assert hashlib.sha256(p.read_bytes()).hexdigest()==h,p
w=read(P/'witness.json'); t=read(pr/'certificate.json'); b=read(db/'witness.json'); o=read(nc/'certificate.json')
for l,h in w['SHA256'].items():
 d,n=l.split('/')
 H((pr if d=='triangle' else db)/n,h)
H(nc/'certificate.json',t['source_sha256'])
Q=np.asarray(o['gram_12'],np.int64)
M=np.asarray(o['sublattice_matrix'],np.int64)
rot=np.eye(12,dtype=np.int64)-M
ns=importlib.util.spec_from_file_location('shells',nc/'verify.py')
lm=importlib.util.module_from_spec(ns);ns.loader.exec_module(lm)
sh=list(lm.enumerate_short(Q.tolist()))+[(tuple([0]*12),0)]
D={}
for v,q in sh:
 lab=sum((((v[i]+v[j])%3)*3**k) for k,(i,j) in enumerate(zip(o['label_first'],o['label_second'])))
 D.setdefault(lab,[]).append((q,v))
assert len(D)==729
N0={k:min(q for q,v in vals) for k,vals in D.items()}
S={k:np.array([v for q,v in D[k] if q==N0[k]],dtype=np.int64) for k in D}
assert Counter((N0[k],len(S[k])) for k in D)=={(0,1):1,(4,3):252,(6,18):224,(8,81):252}
T={}
for k,s in S.items():
 if N0[k]==4:
  T[k]=np.array([s[0],s[0]@rot,s[0]@rot@rot])
  assert set(map(tuple,T[k]))==set(map(tuple,s))
cells=o['outer_cosets']
rt={tuple(c)+p for c in cells if [N0[k] for k in c]==[4]*3 for p in itertools.product(range(3),repeat=3)}
for obj,ra,ad in [(t,'removed_points','added_points'),(b,'remove_rows',None)]:
 rem=set(map(tuple,obj[ra])); assert len(rem)==len(obj[ra]) and rem<=rt; rt-=rem
 if ad:
  add=set(map(tuple,obj[ad])); assert len(add)==len(obj[ad]) and not (add&rt);rt|=add
 else:
  add={tuple(c)+p for c in obj['insert_full_triples'] for p in itertools.product(range(3),repeat=3)}
  assert not add&rt; rt|=add
assert len(rt)==99225
def root(row):
 return np.concatenate([T[k][a] for k,a in zip(row[:3],row[3:])])
ax=np.array([v for v,q in sh if q==4],dtype=np.int64)@M
base=[]; ids=[]; z=np.zeros(12,dtype=np.int64)
for i in range(3):
 for v in ax:
  A=[z]*3;A[i]=v;base.append(np.concatenate(A));ids.append(None)
moved={(a,tuple(b)) for a,b in t['switched_triangles']}
def neg(n):
 return sum((-(n//3**j)%3)*3**j for j in range(6))
for i,c in enumerate(cells[1:],1):
 q=[N0[j] for j in c]
 if q==[4,4,4]: continue
 if sorted(q)==[0,6,6]:
  for a in itertools.product(*(S[k] for k in c)):
   base.append(np.concatenate(a));ids.append(None)
 else:
  assert sorted(q)==[0,4,8]
  a,h,e=q.index(4),q.index(8),q.index(0)
  for j,v in enumerate(S[c[h]]):
   switch=(i,tuple(v)) in moved
   s,tag=(e,neg(c[a])) if switch else (a,c[a])
   for phase,p in enumerate(T[tag]):
    A=[z]*3; A[s]=p;A[h]=v
    base.append(np.concatenate(A)); ids.append((i,j,phase))
assert len(base)==403704 and len(ids)==len(base)
assert len(set(k for k in ids if k is not None))==sum(k is not None for k in ids)
rold=sorted(rt)
old=np.vstack([np.array(base),np.array([root(row) for row in rold])])
assert len(old)==502929 and len(np.unique(old,axis=0))==502929
adds=[]
for a,b,c,mask in w['add']:
 assert 0<mask<2**27
 adds += [(a,b,c,i%3,i//3%3,i//9) for i in range(27) if mask>>i&1]
assert len(adds)==len(set(adds)) and not (rt&set(adds))
new=np.asarray([root(row) for row in adds])
G=np.zeros((36,36),dtype=np.int64)
for i in range(3):
 G[12*i:12*i+12,12*i:12*i+12]=Q
Ob=old@G
assert np.all(np.sum(old*Ob,axis=1)==12)
Nb=new@G
Pnew=Nb@new.T
assert np.all(Pnew.diagonal()==12)
np.fill_diagonal(Pnew,-99)
assert Pnew.max()==6
del_req={(a,b,i) for a,b,mask in w['remove_base'] for i in range(3) if mask>>i&1}
assert len(del_req)==len([1 for a,b,mask in w['remove_base'] for i in range(3) if mask>>i&1])
bad=np.zeros(len(old),dtype=bool)
actual_max=-100
# Floating point only accelerates literal dot products: all inputs are exact
# integers, and any evaluation order is exact by the absolute partial-sum cap.
cap=36*int(np.abs(Nb).max())*int(np.abs(old).max())
assert cap<2**53
Nbf=Nb.astype(np.float64)
for j in range(0,len(old),8192):
 x=old[j:j+8192].astype(np.float64)
 D=Nbf@x.T
 bad[j:j+len(x)]=(D>6).any(axis=0)
 d=D[:,~bad[j:j+len(x)]]
 if d.size: actual_max=max(actual_max,int(d.max()))
actual_ids={ids[i] for i in np.flatnonzero(bad[:len(base)])}
assert actual_ids<=del_req
actual_root_ids={rold[i] for i in np.flatnonzero(bad[len(base):])}
root_req=set(map(tuple,w['remove_root'])); assert actual_root_ids<=root_req<=rt
assert len(root_req)==len(w['remove_root'])
assert del_req<=set(ids)
assert sum(bad)==len(actual_ids)+len(actual_root_ids)
bad[:len(base)] |= np.asarray([i in del_req for i in ids])
bad[len(base):] |= np.asarray([r in root_req for r in rold])
assert sum(bad)==len(del_req)+len(root_req)
print('Exact base blocker subset',len(actual_ids),'with additional planned deletions',len(del_req-actual_ids),flush=True)
assert actual_max==6
# Independently expand the recovered two-block points and compare them
# literally to every retained point, including the whole new root triple bank.
recover=[]
for own,ii,small,label,mask in w['restore_base']:
 assert 0<own<len(cells) and 0<mask<8
 c=cells[own]; typ=[N0[k] for k in c]
 assert sorted(typ)==[0,4,8]
 h=typ.index(8)
 assert 0<=ii<len(S[c[h]]) and small in range(3) and small!=h
 assert all((own,ii,j) in del_req for j in range(3))
 for j in range(3):
  if mask>>j&1:
   A=[z]*3;A[h]=S[c[h]][ii];A[small]=T[label][j]
   recover.append(np.concatenate(A))
Z=np.asarray(recover)
assert len(Z)>0 and len(np.unique(Z,axis=0))==len(Z)
BZ=Z@G
C=BZ@Z.T
assert np.all(C.diagonal()==12)
np.fill_diagonal(C,-99)
assert C.max()<=6
assert (BZ@new.T).max()<=6
cap2=36*int(np.abs(BZ).max())*int(np.abs(old).max())
assert cap2<2**53
checks_restore=0
max_restore=-99
for j in range(0,len(old),8192):
 x=old[j:j+8192];x=x[~bad[j:j+len(x)]].astype(np.float64)
 checks_restore+=len(Z)*len(x)
 if len(x): max_restore=max(max_restore,int((BZ.astype(np.float64)@x.T).max()))
assert max_restore<=6
assert len(old)+len(new)+len(Z)-sum(bad)==w['points']
print('PASS recovered',len(Z),'with max',max_restore,'in',checks_restore,'old products,',len(Z)*len(new),'new root products and',len(Z)*(len(Z)-1)//2,'within recovered pairs; partial sum cap',cap2)
print('PASS independent literal expansion: exact max',actual_max,'over',len(old)*len(new),'new-old products; float64 absolute partial-sum cap',cap,';', int(w['points']), 'with',len(adds),'added;',len(del_req),'base deletions and',len(root_req),'root deletions, actual root blockers',len(actual_root_ids))
