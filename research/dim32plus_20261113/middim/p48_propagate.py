import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np,gzip,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).parent;P=ROOT/'verifications/improved/dim49-60-p48-pointwise-refill'
new=np.load(OUT/'p48_signed_class_lns.npz')['points'].astype(np.int64);q=np.load(P/'heads49.npz');G=q['gram'];baseline=np.r_[q['base'],q['extra_points']];baseline=np.r_[baseline,-baseline]
oldset=set(map(tuple,baseline));newset=set(map(tuple,new));rem=np.array(sorted(oldset-newset),dtype=np.int64);add=np.array(sorted(newset-oldset),dtype=np.int64)
print('template removed',len(rem),'added',len(add),flush=True)
# Canonical line bytes, and sign relative to that canonical representative.
def canonical(X):
 X=np.array(X,dtype=np.int64).reshape(-1,48)
 j=np.argmax(X!=0,axis=1);s=np.sign(X[np.arange(len(X)),j]);return X*s[:,None],s

def keys(X):return np.ascontiguousarray(X,dtype=np.int16).view(np.dtype((np.void,96))).ravel()

def abskeys(X):
 x,s=canonical(X);assert abs(x).max()<32768;return keys(x),s
poles=[2,6,12,24,40,72,126,240,306,510,604,841]
for d in range(49,61):
 st=time.time();a=np.load(P/f'heads{d}.npz');B=a['base'];T=a['transforms'];bits=np.unpackbits(a['retained_bits'],axis=1)[:,:len(B)].astype(bool);E=a['extra_points'];O=a['extra_offsets'];m=len(T)
 oldlocals=[np.r_[B[bits[i]],E[O[i]:O[i+1]]] for i in range(m)]
 oldglobal=np.concatenate([canonical(x@t)[0].astype(np.int16) for x,t in zip(oldlocals,T)])
 kk=keys(oldglobal);sorter=np.argsort(kk);sortedkeys=kk[sorter]
 assert len(np.unique(sortedkeys))==len(sortedkeys)
 owner=np.repeat(np.arange(m),[len(x) for x in oldlocals]);owner=owner[sorter]
 # Per signed mother point state can now differ from its antipode.
 deleted=set(); inserted=set();changes=[];totalgain=0
 for i in range(m):
  old=oldlocals[i];old=np.r_[old,-old]
  lookup={tuple(x):j for j,x in enumerate(old)}
  remove=[x for x in rem if tuple(x) in lookup]
  removing=set(map(tuple,remove));retained=np.array([x for x in old if tuple(x) not in removing])
  # Complete class replay against all retained old vectors (integer products exact).
  available=[]
  ag=add@G
  for j,x in enumerate(add):
   if tuple(x) in lookup:continue
   if np.max(ag[j]@retained.T)>2:continue
   available.append(x)
  if not available:continue
  available=np.array(available)
  global_add=available@T[i];ka,sa=abskeys(global_add)
  keep=[]
  for j,(k,s) in enumerate(zip(ka,sa)):
   key=(k.tobytes(),int(s));pos=np.searchsorted(sortedkeys,k)
   occupied=(pos<len(sortedkeys) and sortedkeys[pos]==k and key not in deleted) or key in inserted
   if not occupied:keep.append(j)
  available=available[keep]
  if len(available)<=len(remove):continue
  # Pairs among added rows are inherited from the newly verified signed template.
  assert np.max((available@G@available.T)-6*np.eye(len(available),dtype=np.int64))<=2
  global_removed=np.array(remove,dtype=np.int64).reshape(-1,48)@T[i];kr,sr=abskeys(global_removed)
  for k,s in zip(kr,sr):deleted.add((k.tobytes(),int(s)))
  ka,sa=abskeys(available@T[i])
  for k,s in zip(ka,sa):inserted.add((k.tobytes(),int(s)))
  changes.append(dict(class_index=i,remove=np.array(remove).reshape(-1,48).tolist(),add=available.tolist()))
  totalgain+=len(available)-len(remove)
 result=dict(dimension=d,bound=52416000+2*len(oldglobal)+totalgain+poles[d-49],baseline=52416000+2*len(oldglobal)+poles[d-49],gain=totalgain,changed_classes=len(changes),signed_removed=sum(len(c['remove']) for c in changes),signed_added=sum(len(c['add']) for c in changes),changes=changes,seconds=time.time()-st)
 (OUT/f'p48_signed_delta{d}.json.gz').write_bytes(gzip.compress(json.dumps(result,separators=(',',':')).encode(),compresslevel=9,mtime=0))
 print({k:v for k,v in result.items() if k!='changes'},flush=True)
