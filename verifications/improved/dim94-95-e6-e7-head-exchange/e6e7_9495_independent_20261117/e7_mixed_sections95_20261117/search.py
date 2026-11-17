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

P38=pool;M38=[mask(v) for v in pool];base38=P38[0];m0=M38[0]
labels32=[];start32=[]
for a,b in itertools.product(range(8),repeat=2):
 if a%4 not in (2,3):continue
 h=np.zeros(8,dtype=np.int16);h[a%4]=2*(-1)**(a//4);h[4+b%4]=2*(-1)**(b//4)
 h[[1,4]]=h[[4,1]];h[4]*=-1
 assert h[0]+h[4]==0
 labels32.append([a,b]);start32.append(lookup[tuple(h)])
pool=[np.array(start32,dtype=np.int16)];seen={mask(start32)};q=0
while q<len(pool):
 old=pool[q];q+=1
 for perm in perms:
  v=perm[old];m=mask(v)
  if m not in seen:seen.add(m);pool.append(v)
P32=pool;M32=[mask(v) for v in pool];print('ORBIT32',len(pool),'bestpairs',max((m0|m).bit_count() for m in M38),max((m0|m).bit_count() for m in M32),flush=True)
best=(0,None)
for i,m in enumerate(M38):
 u=m0|m
 for j,n in enumerate(M32):
  size=(u|n).bit_count()
  if size>best[0]:best=(size,('38+38+32',i,j));print('BEST',best,flush=True)
for i,m in enumerate(M32):
 u=m0|m
 if u.bit_count()+32<=best[0]:continue
 for j,n in enumerate(M32[i+1:],i+1):
  size=(u|n).bit_count()
  if size>best[0]:best=(size,('38+32+32',i,j));print('BEST',best,flush=True)
(P/'result.json').write_text(json.dumps(dict(best=best,orbit38=len(P38),orbit32=len(P32)))+'\n')
print('FINAL',best)
size,(kind,i,j)=best
if kind=='38+38+32':tables=[base38,P38[i],P32[j]];prefixes=[labels,labels,labels32]
else:tables=[base38,P32[i],P32[j]];prefixes=[labels,labels32,labels32]
used=set();selected=[]
for color,(T,L) in enumerate(zip(tables,prefixes)):
 for lab,ri in zip(L,T):
  ri=int(ri)
  if ri in used:continue
  used.add(ri);selected.append(dict(color=color,prefix=lab,head=R[ri].tolist()))
assert len(selected)==size
(P/'candidate.json').write_text(json.dumps(dict(union_size=size,dense_count=size*8**10,source_prefixes_per_color=prefixes,full_rotated_heads=[R[T].tolist() for T in tables],selected=selected),indent=2)+'\n')
