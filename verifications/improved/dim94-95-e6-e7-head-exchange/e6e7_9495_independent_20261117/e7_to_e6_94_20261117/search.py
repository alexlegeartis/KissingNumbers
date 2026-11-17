from pathlib import Path
import numpy as np,itertools,json,collections
P=Path(__file__).resolve().parent
# Reconstruct integer E7 norm8 roots and exact Weyl reflection permutations.
R=[]
for i,j in itertools.combinations(range(8),2):
 for a,b in itertools.product((-2,2),repeat=2):
  v=[0]*8;v[i]=a;v[j]=b
  if v[0]+v[4]==0:R.append(v)
for s in itertools.product((-1,1),repeat=8):
 if sum(x<0 for x in s)%2==0 and s[0]+s[4]==0:R.append(s)
R=np.array(R,dtype=np.int16);assert R.shape==(126,8);lookup={tuple(v):i for i,v in enumerate(R)}
perms=[]
for r in R:
 dot=R@r;assert np.all(dot%4==0);perms.append([lookup[tuple(v)] for v in R-(dot//4)[:,None]*r])
perms=np.array(perms,dtype=np.int16)
labels=[];start=[]
for a,b in itertools.product(range(8),repeat=2):
 h=np.zeros(8,dtype=np.int16);h[a%4]=2*(-1)**(a//4);h[4+b%4]=2*(-1)**(b//4)
 if h[0]+h[4]==0:labels.append([a,b]);start.append(lookup[tuple(h)])
def mask(a):
 m=0
 for i in a:m|=1<<int(i)
 return m
pool=[np.array(start,dtype=np.int16)];seen={mask(start)};q=0
while q<len(pool):
 old=pool[q];q+=1
 for perm in perms:
  v=perm[old];m=mask(v)
  if m not in seen:seen.add(m);pool.append(v)
assert len(pool)==630
normal=np.array([1,2,0,0,-1,0,0,0]);E6mask=mask(np.flatnonzero(R@normal==0));assert E6mask.bit_count()==72
sets=[mask(v)&E6mask for v in pool];hist=dict(collections.Counter(s.bit_count() for s in sets));print('sections',len(sets),hist,flush=True)
order=sorted(range(len(pool)),key=lambda i:sets[i].bit_count(),reverse=True);best=(0,None)
for ii,i in enumerate(order):
 for jj in range(ii+1,len(order)):
  j=order[jj];u=sets[i]|sets[j]
  if u.bit_count()+sets[order[0]].bit_count()<=best[0]:continue
  for k in order[jj+1:]:
   if u.bit_count()+sets[k].bit_count()<=best[0]:continue
   n=(u|sets[k]).bit_count()
   if n>best[0]:best=(n,[i,j,k]);print('BEST',n,flush=True)
n,choice=best;used=set();selected=[]
for c,idx in enumerate(choice):
 for p,ri in zip(labels,pool[idx]):
  ri=int(ri)
  if not(E6mask>>ri&1) or ri in used:continue
  used.add(ri);selected.append(dict(color=c,prefix=p,head=R[ri].tolist()))
obj=dict(union_size=n,dense_count=n*8**10,E6_normal=normal.tolist(),section_histogram=hist,orbit_count=630,source_prefixes=labels,full_rotated_heads=[R[pool[i]].tolist() for i in choice],selected=selected)
(P/'candidate.json').write_text(json.dumps(obj,indent=2)+'\n');print('FINAL',n,n*8**10)
