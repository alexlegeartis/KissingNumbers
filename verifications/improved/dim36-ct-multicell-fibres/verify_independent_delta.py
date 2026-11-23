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
w=read(Path(sys.argv[1])); t=read(pr/'certificate.json'); b=read(db/'witness.json'); o=read(nc/'certificate.json')
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

def assemble(W):
 rd=set(map(tuple,W['remove_root']))
 assert len(rd)==len(W['remove_root']) and rd<=rt
 ba=list((a,b,i) for a,b,mask in W['remove_base'] for i in range(3) if mask>>i&1)
 assert all(0<mask<8 for _,__,mask in W['remove_base'])
 assert len(set(ba))==len(ba) and set(ba)<=set(ids)
 ix={v:i for i,v in enumerate(ids) if v is not None}
 bm=np.ones(len(base),dtype=bool)
 bm[[ix[v] for v in ba]]=False
 ad=[]
 for a,b,c,mask in W['add']:
  assert 0<mask<2**27
  ad.extend((a,b,c,i%3,i//3%3,i//9) for i in range(27) if mask>>i&1)
 assert len(ad)==len(set(ad)) and not rt&set(ad)
 rs=sorted((rt-rd)|set(ad))
 ae=[]
 for a,b,c,d,mask in W['restore_base']:
  assert 0<a<len(cells) and 0<mask<8
  typ=[N0[i] for i in cells[a]]
  assert sorted(typ)==[0,4,8]
  h=typ.index(8)
  assert 0<=b<len(S[cells[a][h]]) and c in range(3) and c!=h
  assert all((a,b,k) in ix and not bm[ix[a,b,k]] for k in range(3))
  ae.extend((a,b,c,d,i) for i in range(3) if mask>>i&1)
 assert len(ae)==len(set(ae))
 assert bm.sum()+len(rs)+len(ae)==W['points']
 return bm,rs,sorted(ae)
def rbase(desc):
 a,b,c,d,i=desc
 g=next(j for j,k in enumerate(cells[a]) if N0[k]==8)
 row=[z]*3;row[g]=S[cells[a][g]][b];row[c]=T[d][i]
 return np.concatenate(row)
v=read(Path(sys.argv[2]))
assert v['SHA256']==w['SHA256']
B0,R0,A0=assemble(w)
B1,R1,A1=assemble(v)
d0=(~B0)[B1]
rr=set(R0);d1=[x not in rr for x in R1]
aa=set(A0);d2=[x not in aa for x in A1]
di=np.flatnonzero(np.r_[d0,d1,d2])
F=np.vstack([old[:len(base)][B1],np.array([root(x) for x in R1]),np.array([rbase(x) for x in A1])])
assert len(F)==v['points']
print("Ready delta by type",sum(d0),sum(d1),sum(d2),"from previously independently known",w['points'],"to",v['points'],flush=True)
G=np.kron(np.eye(3,dtype=np.int64),Q)
B=F[di]@G
assert np.all(np.sum(F[di]*B,axis=1)==12)
assert 36*int(np.abs(B).max())*int(np.abs(F).max())<2**53
B=B.astype(np.float64)
for j in range(0,len(F),4096):
 X=F[j:j+4096].astype(np.float64)
 C=B@X.T
 I=np.flatnonzero((di>=j)&(di<j+len(X)))
 assert np.all(C[I,di[I]-j]==12)
 C[I,di[I]-j]=-99
 assert C.max()<=6
print("PASS exact delta against all rows; compared",len(di)*len(F),"products; count",v['points'])
