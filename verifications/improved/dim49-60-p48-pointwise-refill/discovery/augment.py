"""Regenerate ownership and refills from saved initial isometries and a pool.
Usage: python augment.py OUTPUT_DIR [DIMENSIONS ...]
The final witness files supply only initial transforms/IDs; masks and extras
are rebuilt from scratch, independently of their recorded final contents.
"""
import os
os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["OMP_NUM_THREADS"]="1"
import numpy as np,json,time,sys
from pathlib import Path
PKG=Path(__file__).resolve().parent.parent
if len(sys.argv)<2:raise SystemExit('Usage: python augment.py OUTPUT_DIR [DIMENSIONS ...]')
D=Path(sys.argv[1]).resolve();D.mkdir(parents=True,exist_ok=True)
S=PKG/'source'
z=np.load(S/'family_refilled.npz');B=z['base'].astype(np.int64);G=z['gram'].astype(np.int64);pool_archive=np.load(D/'low_conflict_pool.npz');pool_images=int(pool_archive['images']);P=pool_archive['points'].astype(np.int64)
# Recompute exact bounded integer conflicts; exclude existing base lines.
assert 48*int(abs(P@G).max())*int(abs(B).max())<2**63
Cs=[];keep=[]
for st in range(0,len(P),256):
 assert 48*int(abs(P[st:st+256]@G).max())*int(abs(B).max())<2**53
 H=(P[st:st+256]@G).astype(float)@B.astype(float).T
 for j,h in enumerate(H):
  if np.max(abs(h))>=6:continue
  cc=np.flatnonzero(abs(h)>2);assert len(cc)<=3;Cs.append(cc);keep.append(st+j)
P=P[keep];order=sorted(range(len(P)),key=lambda i:(len(Cs[i]),tuple(P[i])));P=P[order];Cs=[Cs[i] for i in order]
assert np.all(((P@G)*P).sum(1)==6)
ms=[1,3,6,12,20,36,63,120,153,255,302,420];poles=[2,6,12,24,40,72,126,240,306,510,604,841]
def canon(X):
 j=np.argmax(X!=0,axis=1);return X*np.sign(X[np.arange(len(X)),j])[:,None]
def key(x):return x.astype(np.int8).tobytes()
for dim in ([int(x) for x in sys.argv[2:]] or list(range(49,61))):
 if dim not in range(49,61):raise ValueError('Dimensions must be between 49 and 60')
 start=time.monotonic();m=ms[dim-49]
 q=np.load(PKG/f'heads{dim}.npz',allow_pickle=False);ids=q['class_ids'];T=q['transforms']
 masks=np.zeros((m,len(B)),bool);owner={};alt={};X=[]
 for i,t in enumerate(T):
  assert np.array_equal(t@G@t.T,G);x=canon(B@t);assert np.max(abs(x))<128;X.append(x.astype(np.int8))
  for j,row in enumerate(x):
   k=key(row)
   if k not in owner:owner[k]=(i,j);masks[i,j]=1
   else:
    if k not in alt:alt[k]=[owner[k]]
    alt[k].append((i,j))
 old=len(owner);extras=[[] for _ in range(m)];extraGs=[None]*m;moves=0
 # Three passes allow later ownership changes to unlock previously blocked candidates.
 for sweep in range(3):
  gains=0
  for i in range(m):
   C=canon(P@T[i]);assert np.max(abs(C))<128
   for a,p in enumerate(P):
    k=key(C[a])
    if k in owner:continue
    if extras[i] and np.any(abs(extraGs[i]@p)>2):continue
    blocked=[int(j) for j in Cs[a] if masks[i,j]];plan=[];ok=True
    for j in blocked:
     found=None
     for ni,nj in alt.get(X[i][j].tobytes(),[]):
      if ni==i:continue
      if extras[ni] and np.any(abs(extraGs[ni]@B[nj])>2):continue
      found=(ni,nj);break
     if found is None:ok=False;break
     plan.append((j,*found))
    if not ok:continue
    for j,ni,nj in plan:
     assert masks[i,j] and not masks[ni,nj];masks[i,j]=False;masks[ni,nj]=True;owner[X[i][j].tobytes()]=(ni,nj);moves+=1
    extras[i].append(p);extraGs[i]=np.array(extras[i])@G;owner[k]=(i,-1);gains+=1
   if (i+1)%50==0:print(dim,'sweep',sweep,'class',i+1,'gain',len(owner)-old,'moves',moves,'seconds',round(time.monotonic()-start,1),flush=True)
  print(dim,'sweep_done',sweep,'new',gains,flush=True)
  if not gains:break
 offsets=np.r_[0,np.cumsum([len(e) for e in extras])];E=np.concatenate([np.array(e,dtype=np.int64).reshape(-1,48) for e in extras]);total=int(masks.sum()+len(E));assert total==len(owner)
 out=D/f'heads{dim}.npz';np.savez_compressed(out,class_ids=ids,transforms=T,retained_bits=np.packbits(masks,axis=1),extra_points=E,extra_offsets=offsets,base=B,gram=G)
 result=dict(dimension=dim,selected_lines=total,baseline_lines=old,added_lines=total-old,moves=moves,bound=52416000+2*total+poles[dim-49],poles=poles[dim-49],classes=m,pool_images=pool_images,pool_candidates=len(P),seconds=time.monotonic()-start)
 (D/f'result{dim}.json').write_text(json.dumps(result,indent=2));print('RESULT',result,flush=True)
 del owner,alt,X
