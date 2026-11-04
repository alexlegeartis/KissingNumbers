"""Seeded competitive DGLS ring-coset and coordinate-dependent binary-mask screen.
All rejection witnesses are actual parent messages; no zero-overlap rejection.
"""
from pathlib import Path
from itertools import combinations
from collections import Counter,deque
import json,hashlib,time
import numpy as np
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2]
RAW=ROOT/'verifications/improved/dim88-96-z8-golay/data/z8_generator.json'
SEED=20261104096
T=np.array([sum((((a>>2)&1)^((a&1)*(x&1))^(((a>>1)&1)*((x>>1)&1)))<<x for x in range(4))for a in range(8)],dtype=np.uint8)
PC=np.array([i.bit_count()for i in range(16)],dtype=np.uint8)
H=np.array([0,2,2,2,4,2,2,2],dtype=np.int16)

def info_set(G,permutation):
 A=G.copy().astype(np.int64);R=np.eye(12,dtype=np.int64);pivots=[];at=0
 for col in permutation:
  candidates=np.flatnonzero(A[at:,col]%2)
  if not len(candidates):continue
  j=at+int(candidates[0]);A[[at,j]]=A[[j,at]];R[[at,j]]=R[[j,at]]
  inverse=pow(int(A[at,col]),-1,8);A[at]=A[at]*inverse%8;R[at]=R[at]*inverse%8
  for i in range(12):
   if i!=at:
    c=int(A[i,col]);A[i]=(A[i]-c*A[at])%8;R[i]=(R[i]-c*R[at])%8
  pivots.append(int(col));at+=1
  if at==12:break
 if at!=12 or not np.array_equal(R@G%8,A):raise ValueError('invalid information set')
 return np.array(pivots),A,R

def canonical(V,base):
 I,A,R=base
 return (V-V[:,I]@A)%8

def ring_screen(V,infos):
 best=np.full(len(V),97,dtype=np.int16);messages=np.zeros((len(V),12),dtype=np.int16)
 for I,A,R in infos:
  residual=(V-V[:,I]@A)%8;d=H[residual].sum(axis=1)
  improve=d<best;best[improve]=d[improve];messages[improve]=(V[improve][:,I]@R)%8
 return best,messages

def binary_screen(M,infos,rng,rounds=1):
 best=np.full(len(M),97,dtype=np.int16);messages=np.zeros((len(M),12),dtype=np.int16);which=np.zeros(len(M),dtype=np.int16)
 # Reject by nearest-to-either-existing-layer witnesses only. A survivor remains
 # unverified because full crossdistance also varies the second parent word.
 local=PC[T[None,:]^np.arange(16,dtype=np.uint8)[:,None]]
 nearest=[np.flatnonzero(local[m]==local[m].min())for m in range(16)]
 for layer in (0,1):
  target=M^(1 if layer else 0)
  for I,A,R in infos:
   for repeat in range(rounds):
    chosen=np.zeros((len(M),12),dtype=np.int64)
    for j,col in enumerate(I):
     for symbol in range(16):
      ids=np.flatnonzero(target[:,col]==symbol)
      if len(ids):chosen[ids,j]=rng.choice(nearest[symbol],len(ids))
    W=chosen@A%8;d=PC[T[W]^target].sum(axis=1)
    improve=(d<best)&(d>0);best[improve]=d[improve];messages[improve]=(chosen[improve]@R)%8;which[improve]=layer
 return best,messages,which

def binary_covering(G,base):
 I,A,R=base;outside=[j for j in range(24)if j not in I]
 columns=[]
 for j in range(24):
  e=np.zeros((1,24),dtype=np.int64);e[0,j]=1
  s=canonical(e,base)[0,outside]%2
  columns.append(sum(int(v)<<k for k,v in enumerate(s)))
 distances=[-1]*4096;distances[0]=0;Q=deque([0])
 while Q:
  a=Q.popleft()
  for b in columns:
   if distances[a^b]<0:distances[a^b]=distances[a]+1;Q.append(a^b)
 return outside,distances

def main():
 start=time.monotonic();rng=np.random.default_rng(SEED);G=np.array(json.loads(RAW.read_text())['generator'],dtype=np.int64)
 base=info_set(G,list(range(24)));infos=[base]+[info_set(G,rng.permutation(24))for _ in range(127)]
 outside,cover=binary_covering(G,base);deep=np.flatnonzero(np.array(cover)==max(cover))
 ring=[];saved_ring_offsets=[];saved_ring_messages=[];saved_ring_distances=[]
 # Uniform canonical ring cosets; deep binary-residue cosets; quadratic block offsets.
 for kind in ('uniform_ring_syndrome','deep_residue_lift','quadratic_coordinate_offsets'):
  if kind=='uniform_ring_syndrome':
   V=np.zeros((4096,24),dtype=np.int64);V[:,outside]=rng.integers(0,8,(len(V),12))
  elif kind=='deep_residue_lift':
   V=np.zeros((4096,24),dtype=np.int64);labels=rng.choice(deep,len(V));bits=(labels[:,None]>>np.arange(12))&1
   V[:,outside]=bits+2*rng.integers(0,4,(len(V),12))
  else:
   labels=np.arange(24);bits=((labels[:,None]>>np.arange(5))&1).T
   features=np.vstack([np.ones(24,dtype=int),bits,*[bits[i]*bits[j]for i in range(5)for j in range(i)]])
   V=canonical(rng.integers(0,8,(4096,len(features)))@features%8,base)
  V=np.unique(V,axis=0);V=V[np.any(V,axis=1)]
  best,coeff=ring_screen(V,infos[:32]);survive=best>=24
  saved_ring_offsets.append(V);saved_ring_messages.append(coeff);saved_ring_distances.append(best)
  ring.append(dict(family=kind,candidates=len(V),rejected=int((best<24).sum()),survivors=int(survive.sum()),minimum_witness_distance_histogram=dict(Counter(map(int,best))),best_proxy=int(best.max())))
  print('RING',ring[-1],flush=True)
  # Keep actual maximizers for bounded local syndrome ascent and exact witness audit.
  ids=np.argsort(best)[-32:]
  if kind=='uniform_ring_syndrome':elite=V[ids]
  else:elite=np.vstack([elite,V[ids]])
 # Maximin ascent preserves the binary-residue class via even syndrome changes.
 history=[];evaluated=0
 for generation in range(24):
  parents=elite[rng.integers(len(elite),size=512)].copy()
  columns=rng.choice(outside,512);parents[np.arange(512),columns]=(parents[np.arange(512),columns]+rng.choice([2,4,6],512))%8
  V=np.unique(np.vstack([elite,parents]),axis=0);d,c=ring_screen(V,infos[:32]);evaluated+=len(V)
  saved_ring_offsets.append(V);saved_ring_messages.append(c);saved_ring_distances.append(d)
  ids=np.argsort(d)[-64:];elite=V[ids]
  history.append(dict(generation=generation,candidates=len(V),best_proxy=int(d.max()),survivors=int((d>=24).sum())))
 d,c=ring_screen(elite,infos);audit=[]
 for v,distance,msg in zip(elite,d,c):
  audit.append(dict(offset=v.tolist(),parent_message=msg.tolist(),witness_distance=int(distance)))
 print('ASCENT',history[-1],'expanded128best',int(d.max()),'survivors',int((d>=24).sum()),flush=True)
 # Genuine binary translates with nonconstant blockparity, hence disjoint from
 # both existing even-all-block and odd-all-block layers by parity alone.
 binary=[];binary_audit=[];saved_binary_masks=[];saved_binary_messages=[];saved_binary_distances=[];saved_binary_layers=[]
 for family in ('balanced_random_masks','quadratic_block_masks','golay_residue_parity_masks'):
  count=2048
  if family=='balanced_random_masks':
   M=rng.integers(0,16,(count,24),dtype=np.uint8)
  elif family=='quadratic_block_masks':
   labels=np.arange(24);bits=((labels[:,None]>>np.arange(5))&1).T
   F=np.vstack([np.ones(24,dtype=int),bits,*[bits[i]*bits[j]for i in range(5)for j in range(i)]])
   M=np.zeros((count,24),dtype=np.uint8)
   for bit in range(4):M|=(((rng.integers(0,2,(count,len(F)))@F)%2).astype(np.uint8)<<bit)
  else:
   coeff=rng.integers(0,2,(count,12));residue=coeff@G%2
   # Bias parity toward nonconstant low-weight Golay supports, with random local masks.
   M=rng.integers(0,16,(count,24),dtype=np.uint8)
   parity=PC[M]%2;M^=(parity^residue).astype(np.uint8)
  parity=(PC[M]%2).sum(axis=1);M=M[(parity>0)&(parity<24)]
  best,msg,layer=binary_screen(M,infos[:32],rng)
  saved_binary_masks.append(M);saved_binary_messages.append(msg);saved_binary_distances.append(best);saved_binary_layers.append(layer)
  binary.append(dict(family=family,candidates=len(M),rejected=int((best<24).sum()),survivors=int((best>=24).sum()),witness_distance_histogram=dict(Counter(map(int,best))),best_proxy=int(best.max())))
  print('BINARY',binary[-1],flush=True)
  for index in np.argsort(best)[-16:]:binary_audit.append(dict(mask_blocks=M[index].tolist(),parent_message=msg[index].tolist(),existing_layer=int(layer[index]),witness_distance=int(best[index])))
 np.savez_compressed(OUT/'witnesses.npz',ring_offsets=np.vstack(saved_ring_offsets).astype(np.uint8),ring_messages=np.vstack(saved_ring_messages).astype(np.uint8),ring_distances=np.concatenate(saved_ring_distances),binary_masks=np.vstack(saved_binary_masks),binary_messages=np.vstack(saved_binary_messages).astype(np.uint8),binary_distances=np.concatenate(saved_binary_distances),binary_layers=np.concatenate(saved_binary_layers))
 receipt=dict(status='bounded seeded competitive search; survivors are unverified',seed=SEED,generator_sha256=hashlib.sha256(RAW.read_bytes()).hexdigest(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  objective=dict(parent_cardinality=2**36,existing_two_layer_cardinality=2**37,third_disjoint_layer_cardinality=3*2**36,repository_comparison=142162362560,required_cross_distance=24),
  ring_coset_feasibility='Gray(C+v) has parent internal distance by isometry; disjoint even-block ring cosets plus existing odd-all-block translate. A third ring coset needs exact min homogeneous distance>=24.',
  binary_coset_feasibility='C xor mask has internal distance; nonconstant blockparity proves disjointness from both existing layers. Nearest-point witnesses reject; survivors do not certify full crossdistance.',
  residue_covering_radius=max(cover),residue_covering_histogram=dict(Counter(cover)),information_sets=len(infos),ring_results=ring,ascent_history=history,ascent_candidates=evaluated,ring_elite_witnesses=audit,binary_results=binary,binary_elite_witnesses=binary_audit,elapsed_seconds=time.monotonic()-start)
 (OUT/'results.json').write_text(json.dumps(receipt,indent=2)+'\n');print('elapsed',receipt['elapsed_seconds'],flush=True)
if __name__=='__main__':main()
