"""Exact 39D QR pentapolar certificate replay using Python integers."""
from itertools import combinations,product
import gzip, hashlib, json, math
from pathlib import Path
P=Path(__file__).resolve().parent
C=json.loads((P/'certificate.json').read_text())
D=json.loads((P/'parent39.json').read_text())
B=json.loads((P/'base_certificate.json').read_text())
marks=C['marked_axes']; un=[i for i in range(39) if i not in marks]; assert len(marks)==2 and marks==sorted(set(marks)) and len(un)==37
L=C['scale']; assert type(L)==int and L>0
hh=C['heads']; assert len(hh)==5 and all(len(h)==2 and all(type(x)==int for x in h) for h in hh)
g=[int(s,16) for s in D['generators_hex']]; orig=[int(s,16) for s in D['representatives_hex']]
R={}
for x in g:
 while x:
  i=x.bit_length()-1
  if i not in R:R[i]=x;break
  x^=R[i]
assert len(R)==15
free=sorted(set(range(39))-R.keys())
def q(x):
 for i,v in sorted(R.items(),reverse=True):
  if x>>i&1:x^=v
 return sum(((x>>i)&1)<<j for j,i in enumerate(free))
cols=[q(1<<i) for i in un]; Q={}
for j,v in enumerate(cols):
 u=1<<j
 while v:
  i=v.bit_length()-1
  if i not in Q:Q[i]=(v,u);break
  d,h=Q[i];v^=d;u^=h
assert len(Q)==24
def rep(x):
 u=0
 for i,(v,t) in sorted(Q.items(),reverse=True):
  if x>>i&1:x^=v;u^=t
 assert x==0
 return u
def compress(w):
 return sum(((w>>i)&1)<<j for j,i in enumerate(un))
H=[0]
for v in g:H.extend([w^v for w in H])
K=[compress(w) for w in H if not any(w>>i&1 for i in marks)]
assert len(K)==len(set(K))==8192
mdK=min(w.bit_count() for w in K if w)
A=[q(w) for w in orig]; assert len(A)==len(set(A))==10
T=[0,q(1<<marks[0]),q(1<<marks[1]),q((1<<marks[0])|(1<<marks[1]))]
offs=[[x^z for x in A] for z in T]+[[C['extra_syndrome']]]
F=[(i,rep(s)) for i,l in enumerate(offs) for s in l]
assert len(F)==41 and len(set(F))==41
N=[148*L*L+sum(t*t for t in h) for h in hh]
def angular(dot,n,m):
 return dot<=0 or 4*dot*dot<=n*m
def dense(i,j,d):
 v=(148-8*d)*L*L+sum(a*b for a,b in zip(hh[i],hh[j]))
 assert angular(v,N[i],N[j]), ('dense',i,j,d,v)
for i in range(5):dense(i,i,mdK)
tiny=[[37]*5 for _ in range(5)]
for (i,a),(j,b) in combinations(F,2):
 d=min((k^a^b).bit_count() for k in K);dense(i,j,d);tiny[i][j]=min(tiny[i][j],d)
# Retained sparse points are a subset of the 769852 construction's sparse layer.
raw=gzip.decompress((P/'supports_3424.txt.gz').read_bytes())
assert hashlib.sha256(raw).hexdigest()==B['source_sha256']['supports_3424.txt']
rows=raw.decode('ascii').splitlines()
assert len(rows)==3424 and all(len(row)==39 and set(row)<=set('01') and row.count('1')==8 for row in rows)
supports=[[i for i,b in enumerate(row) if b=='1'] for row in rows]
for v in B['support_replacements']:
 assert supports[v['remove_index']]==v['remove'];supports[v['remove_index']]=v['add']
for x in C.get('weighted_support_exchange', []):
 i=x['index']; t=x['replacement']
 assert set(marks)<=set(supports[i]) and not set(marks)<=set(t)
 supports[i]=t
assert all(len(s)==len(set(s))==8 and all(0<=x<39 for x in s) for s in supports)
covered=set()
for s in supports:
 for t in combinations(sorted(s),5):
  assert t not in covered
  covered.add(t)
def bound(i,a):
 val=sum(2*L*abs(v) for j,v in a.items() if j not in marks)
 val+=sum(hh[i][k]*a.get(j,0) for k,j in enumerate(marks))
 assert angular(val,N[i],32),('sparse',i,a,val)
drop=0; octad_points=0
for s in supports:
 proper=[j for j in marks if j in s]
 for signs in product((2,-2),repeat=8):
  if sum(v<0 for v in signs)%2:continue
  a=dict(zip(s,signs))
  if all(j in a and a[j]>0 for j in marks):
   drop+=1;continue
  for i in range(5):bound(i,a)
  octad_points+=1
assert octad_points+drop==3424*128
removed=set();stars=[]
for v in B['heavy_stars']:
 a=v['center']; b=v['leaves']; assert a not in b and len(b)==len(set(b))==4
 for j in b:removed.add(tuple(sorted((a,j))))
 for s in (-4,4):
  for f in product((-2,2),repeat=4):
   if sum(v<0 for v in f)%2:continue
   stars.append({a:s,**dict(zip(b,f))})
assert len(stars)==len({tuple(sorted(a.items())) for a in stars})==64
def dot(a,b):return sum(v*b.get(j,0) for j,v in a.items())
for a,b in combinations(stars,2):assert dot(a,b)<=16
for a in stars:
 assert sum(v*v for v in a.values())==32
 assert all(2*sum(abs(a.get(j,0)) for j in s)<=16 for s in supports)
 for i in range(5):bound(i,a)
roots=0
for a,b in combinations(range(39),2):
 if (a,b) in removed:continue
 for s,t in product((-4,4),repeat=2):
  v={a:s,b:t}
  assert all(dot(v,q)<=16 for q in stars)
  for i in range(5):bound(i,v)
  roots+=1
# Root/octad and root/root for distinct points are automatic from their shape.
total=len(F)*len(K)+octad_points+roots+len(stars)
assert total==C['claimed_size']
print(json.dumps(dict(status='PASS',dimension=39,points=total,improvement=total-769852,dense=len(F)*len(K),octad_points=octad_points,removed_octad_points=drop,roots=roots,stars=len(stars),kernel_minimum_distance=mdK,fibres=len(F),fibres_by_head=list(map(len,offs)),least_intertier_distances=tiny)))
