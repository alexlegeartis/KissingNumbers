"""Deterministic candidate-pool regeneration. Usage: python sieve.py OUTPUT_DIR [IMAGES] [--resume]."""
import os,sys
os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["OMP_NUM_THREADS"]="1"
import numpy as np,time,json,sys
from pathlib import Path
import argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('output_dir',type=Path)
parser.add_argument('images',type=int,nargs='?',default=3500)
parser.add_argument('--resume',action='store_true')
args=parser.parse_args()
if args.images<1:parser.error('images must be positive')
PKG=Path(__file__).resolve().parent.parent
D=args.output_dir.resolve();D.mkdir(parents=True,exist_ok=True)
S=PKG/'source'
z=np.load(S/'family_refilled.npz');B=z['base'].astype(np.int64);G=z['gram'];gs=np.load(S/'p48p_gens.npy').astype(np.int64)
rng=np.random.default_rng(20261108);BG=(B@G).astype(np.float64);Bf=B.astype(float);pool=[];words=[];start=time.monotonic();limit=3
total=args.images
resume=0
if args.resume:
 old=np.load(D/'low_conflict_pool.npz');resume=int(old['images']);assert int(old['seed'])==20261108 and int(old['max_conflicts'])==limit;assert resume<=total;pool=[old['points']];rng.integers(len(gs),size=(resume,30))
for i in range(resume,total):
 word=rng.integers(len(gs),size=30);T=np.eye(48,dtype=np.int64)
 for g in word:
  assert 48*int(abs(T).max())*int(abs(gs[g]).max())<2**63
  T=T@gs[g]
 assert 48*int(abs(B).max())*int(abs(T).max())<2**63
 X=B@T;j=np.argmax(X!=0,axis=1);X*=np.sign(X[np.arange(len(X)),j])[:,None];Y=X.astype(float);ix=np.arange(len(X));bad=np.zeros(len(X),np.int16)
 for st in range(0,len(B),512):
  assert 48*float(abs(Y).max())*float(abs(BG).max())<2**53
  H=Y@BG[st:st+512].T;bad+=(np.abs(H)>2).sum(1).astype(np.int16);ok=bad<=limit;ix=ix[ok];Y=Y[ok];bad=bad[ok]
  if not len(ix):break
 if len(ix):pool.append(X[ix].astype(np.int16))
 if (i+1)%25==0 or i+1==total:
  P=np.unique(np.concatenate(pool),axis=0) if pool else np.empty((0,48),np.int16)
  pool=[P];np.savez_compressed(D/'low_conflict_pool.npz',points=P,images=i+1,seed=20261108,max_conflicts=limit)
  H=P.astype(float)@BG.T;cnt=(abs(H)>2).sum(1);nonself=np.max(abs(H),axis=1)<6;print(i+1,len(P),'novel',int(nonself.sum()),'hist',np.unique(cnt[nonself],return_counts=True),'secs',round(time.monotonic()-start,1),flush=True)
