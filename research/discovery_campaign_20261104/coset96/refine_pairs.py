"""Vary both DGLS parent words for surviving nonconstant-mask candidates."""
from pathlib import Path
from collections import Counter
import json,time,hashlib
import numpy as np
from search import info_set,T,PC,SEED,RAW
OUT=Path(__file__).resolve().parent
start=time.monotonic();rng=np.random.default_rng(SEED+1)
G=np.array(json.loads(RAW.read_text())['generator'],dtype=np.int64)
data=np.load(OUT/'witnesses.npz');ids=np.flatnonzero(data['binary_distances']>=24)
M=data['binary_masks'][ids].copy();best=data['binary_distances'][ids].copy()
first=data['binary_messages'][ids].astype(np.int64);second=np.zeros((len(M),12),dtype=np.int64);layers=data['binary_layers'][ids].copy()
local=PC[T[None,:]^np.arange(16,dtype=np.uint8)[:,None]];nearest=[np.flatnonzero(local[m]==local[m].min())for m in range(16)]
trials=np.zeros(len(M),dtype=np.int64);history=[]
for trial in range(1024):
 active=np.flatnonzero(best>=24)
 if not len(active):break
 I,A,R=info_set(G,rng.permutation(24));layer=trial%2
 other=rng.integers(0,8,(len(active),12));otherwords=other@G%8
 target=M[active]^T[otherwords]^(1 if layer else 0)
 chosen=np.zeros((len(active),12),dtype=np.int64)
 for j,col in enumerate(I):
  for symbol in range(16):
   where=np.flatnonzero(target[:,col]==symbol)
   if len(where):chosen[where,j]=rng.choice(nearest[symbol],len(where))
 words=chosen@A%8;distance=PC[T[words]^target].sum(axis=1)
 better=(distance<best[active])&(distance>0);improved=active[better]
 best[improved]=distance[better];first[improved]=chosen[better]@R%8;second[improved]=other[better];layers[improved]=layer
 trials[active]+=1
 if trial%16==0 or not np.any(best>=24):
  row=dict(trial=trial,remaining=int((best>=24).sum()),best_proxy=int(best.max()));history.append(row);print(row,flush=True)
np.savez_compressed(OUT/'pair_witnesses.npz',source_indices=ids,masks=M,first_messages=first.astype(np.uint8),second_messages=second.astype(np.uint8),layers=layers,distances=best,trials=trials)
receipt=dict(status='bounded exact witness screen; surviving masks remain uncertified',seed=SEED+1,input_candidates=len(M),rejected=int((best<24).sum()),survivors=int((best>=24).sum()),history=history,distance_histogram=dict(Counter(map(int,best))),pair_trials=int(trials.sum()),max_trials_per_candidate=int(trials.max())if len(trials) else 0,elapsed_seconds=time.monotonic()-start,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(OUT/'pair_results.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
