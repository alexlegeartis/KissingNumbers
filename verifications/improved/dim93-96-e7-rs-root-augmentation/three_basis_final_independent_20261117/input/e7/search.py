from pathlib import Path
import numpy as np,itertools,json,random
P=Path(__file__).resolve().parent;roots=[]
for i,j in itertools.combinations(range(8),2):
 for a,b in itertools.product((-2,2),repeat=2):
  v=[0]*8;v[i]=a;v[j]=b
  if v[0]+v[4]==0:roots.append(v)
for s in itertools.product((-1,1),repeat=8):
 if sum(x<0 for x in s)%2==0 and s[0]+s[4]==0:roots.append(s)
R=np.array(roots,dtype=np.int16);assert R.shape==(126,8);lookup={tuple(v):i for i,v in enumerate(R)}
perms=[]
for r in R:
 dot=R@r;assert np.all(dot%4==0);T=R-(dot//4)[:,None]*r;perms.append([lookup[tuple(v)] for v in T])
perms=np.array(perms)
labels=[];idx=[]
for a,b in itertools.product(range(8),repeat=2):
 h=np.zeros(8,dtype=np.int16);h[a%4]=2*(1 if a<4 else -1);h[4+b%4]=2*(1 if b<4 else -1)
 if h[0]+h[4]==0:labels.append([a,b]);idx.append(lookup[tuple(h)])
assert len(idx)==38
rng=random.Random(1441);base=np.array(idx);cur=base.copy();pool=[];seen=set()
def mask(a):
 out=0
 for i in a:out|=1<<int(i)
 return out
m0=mask(base)
for t in range(20000):
 for _ in range(4):cur=perms[rng.randrange(126),cur]
 m=mask(cur)
 if m not in seen:seen.add(m);pool.append((m,cur.copy()))
top=sorted(range(len(pool)),key=lambda i:(m0|pool[i][0]).bit_count(),reverse=True)[:200]
best=(0,None,None)
for i in top:
 u=m0|pool[i][0]
 for j,(m,v) in enumerate(pool):
  n=(u|m).bit_count()
  if n>best[0]:best=(n,i,j);print('BEST',n,'pool',len(pool),flush=True)
n,i,j=best
sets=[base,pool[i][1],pool[j][1]];used=set();selected=[]
for c,ids in enumerate(sets):
 for label,ind in zip(labels,ids):
  if int(ind) in used:continue
  used.add(int(ind));selected.append(dict(color=c,prefix=label,head=R[ind].tolist()))
assert len(selected)==n
obj=dict(union_size=n,dense_count=n*8**10,full_rotated_heads=[R[a].tolist() for a in sets],source_prefixes=labels,selected=selected,pool_size=len(pool),status='finite candidate; verify exact Gram preservation independently')
(P/'candidate.json').write_text(json.dumps(obj,indent=2)+'\n');print('FINAL',n,n*8**10)
