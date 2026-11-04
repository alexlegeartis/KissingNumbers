"""Competitive three-half-coset assembly search, not a construction claim.
D=kernel of a nonzero binary functional on Z8 message parity. Its size is2^35.
Three distance24 D-cosets plus the existing odd-block full C would give2.5*C.
"""
from pathlib import Path
from itertools import combinations
from collections import Counter
import json,time,hashlib
import numpy as np
from search import info_set,canonical,H,SEED,RAW
OUT=Path(__file__).resolve().parent

def restricted_screen(V,G,ell,infos,rng):
 best=np.full(len(V),97,dtype=np.int16);messages=np.zeros((len(V),12),dtype=np.int16)
 for I,A,R in infos:
  chosen=V[:,I].copy();coeff=chosen@R%8;bad=(coeff@ell)%2
  # Toggle an information-coordinate value whose original-message parity
  # functional is odd. Every tried message then lies in D exactly.
  valid=np.flatnonzero((R@ell)%2);j=int(rng.choice(valid))
  for change in (1,3,5,7):
   modified=chosen.copy();modified[:,j]=(modified[:,j]+change*bad)%8
   c=modified@R%8
   if np.any(c@ell%2):raise ValueError('message outside subcode')
   residual=(V-modified@A)%8;d=H[residual].sum(axis=1)
   improve=(d<best)&(d>0);best[improve]=d[improve];messages[improve]=c[improve]
 return best,messages

def main():
 start=time.monotonic();rng=np.random.default_rng(SEED+2);G=np.array(json.loads(RAW.read_text())['generator'],dtype=np.int64)
 base=info_set(G,range(24));outside=[j for j in range(24)if j not in base[0]]
 infos=[info_set(G,rng.permutation(24))for _ in range(32)]
 functionals=[G[:,j]%2 for j in range(24)]+[rng.integers(0,2,12)for _ in range(24)]
 seen=set();rows=[];survivors_saved=[];total=0;pair_total=0;triples=[]
 for ell in functionals:
  key=tuple(map(int,ell))
  if key in seen or not any(key):continue
  seen.add(key);oddrow=int(np.flatnonzero(ell)[0]);e=G[oddrow]
  random=np.zeros((512,24),dtype=np.int64);random[:,outside]=rng.integers(0,8,(512,12))
  even=np.zeros((256,24),dtype=np.int64);even[:,outside]=2*rng.integers(0,4,(256,12))
  torsion=np.zeros((128,24),dtype=np.int64);torsion[:,outside]=4*rng.integers(0,2,(128,12))
  switching=canonical(4*np.eye(24,dtype=np.int64),base)
  offsets=np.unique(np.vstack([random,even,torsion,switching,np.zeros((1,24),dtype=int)]),axis=0)
  # Both D cosets inside every C coset. Exclude D itself; C's odd half is a known control.
  V=np.unique(np.vstack([offsets, (offsets+e)%8]),axis=0);V=V[np.any(V,axis=1)]
  best,msg=restricted_screen(V,G,ell,infos[:12],rng);total+=len(V)
  ids=np.flatnonzero(best>=24);S=V[ids]
  # Keep all surviving candidates, up to128, prioritizing larger proxy distances.
  if len(S)>128:S=V[np.argsort(best)[-128:]]
  pairs=list(combinations(range(len(S)),2));compatible=[]
  if pairs:
   differences=np.array([(S[i]-S[j])%8 for i,j in pairs],dtype=np.int64)
   pd,pm=restricted_screen(differences,G,ell,infos,rng);pair_total+=len(pairs)
   compatible=[pair for pair,d in zip(pairs,pd)if d>=24]
   for i,j in compatible:triples.append(dict(functional=list(key),offsets=[S[i].tolist(),S[j].tolist()],status='unverified triple; sampled witnesses only'))
  row=dict(functional=list(key),offsets_tested=len(V),surviving_against_D=len(ids),pair_differences_tested=len(pairs),unverified_three_coset_candidates=len(compatible),best_proxy=int(best.max()))
  rows.append(row);print(row,flush=True)
  for index in ids[:8]:survivors_saved.append(dict(functional=list(key),offset=V[index].tolist(),proxy=int(best[index])))
 receipt=dict(status='bounded seeded candidate generation; no exact distance certificate for survivors',seed=SEED+2,forecast=dict(half_coset_size=2**35,three_even_half_cosets_plus_odd_parent=5*2**35,comparison=142162362560),functionals=len(rows),offsets_tested=total,pair_differences_tested=pair_total,rows=rows,surviving_offsets_examples=survivors_saved,unverified_triples=triples,elapsed_seconds=time.monotonic()-start,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 (OUT/'half_cosets.json').write_text(json.dumps(receipt,indent=2)+'\n');print('SUMMARY',len(rows),total,pair_total,len(triples),receipt['elapsed_seconds'],flush=True)
if __name__=='__main__':main()
