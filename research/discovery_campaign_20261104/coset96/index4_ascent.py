"""Maximin candidate generation outside the provably insufficient4-torsion family."""
from pathlib import Path
from collections import Counter
import json,time,hashlib
import numpy as np
from search import info_set,canonical,SEED,RAW
from index4 import plans_for,restricted,quotient_dedup,find_clique
OUT=Path(__file__).resolve().parent
start=time.monotonic();rng=np.random.default_rng(SEED+5);G=np.array(json.loads(RAW.read_text())['generator'],dtype=np.int64)
base=info_set(G,range(24));outside=[j for j in range(24)if j not in base[0]]
infos=[info_set(G,rng.permutation(24))for _ in range(256)]
source=json.loads((OUT/'index4_results.json').read_text());families=[source['rows'][i]['functionals']for i in (0,1,8,16,24,31)]
rows=[];frontier=[];badwitnesses=[];candidate_total=0;all_offsets={}
for fi,raw in enumerate(families):
 L=np.array(raw,dtype=np.int64);plans=plans_for(L,infos,rng);enc=L[0]+2*L[1];p,q=next((p,q)for p in range(12)for q in range(p)if enc[p]and enc[q]and enc[p]!=enc[q]);E=[G[p],G[q]]
 V=np.zeros((1024,24),dtype=np.int64);V[:,outside]=rng.integers(0,8,(1024,12));V[512:,outside]=2*rng.integers(0,4,(512,12))
 for e in E:V=(V+rng.integers(0,2,(len(V),1))*e)%8
 elites=V;pool=[];history=[];total=0
 for generation in range(48):
  if generation:
   parents=elites[rng.integers(len(elites),size=512)].copy()
   cols=rng.choice(outside,len(parents));parents[np.arange(len(parents)),cols]=(parents[np.arange(len(parents)),cols]+rng.integers(1,8,len(parents)))%8
   if generation%3==0:parents=(parents+rng.integers(0,2,(len(parents),1))*E[generation%2])%8
   V=quotient_dedup(np.vstack([elites,parents]),L,base)
  residual=canonical(V,base);V=V[np.any(residual%4,axis=1)]
  # Rotate some query sets to reduce optimization against one frozen oracle.
  ids=list(range(4))+list(range(4+(generation//8)*4,8+(generation//8)*4))
  d,c=restricted(V,L,[plans[i]for i in ids]);total+=len(V)
  if np.any(d>=24):pool.append(V[d>=24])
  eliteids=np.argsort(d)[-64:];elites=V[eliteids]
  history.append(dict(generation=generation,offsets=len(V),best_minimum_upper_estimate=int(d.max()),upper_estimate_at_least24=int((d>=24).sum())))
 pool=quotient_dedup(np.vstack(pool+[elites]),L,base)
 d,c=restricted(pool,L,plans[32:160]);survivor=np.flatnonzero(d>=24);S=pool[survivor];rejected=np.flatnonzero(d<24)
 for i in rejected[:32]:badwitnesses.append(dict(functionals=L.tolist(),offset_difference=pool[i].tolist(),parent_message=c[i].tolist(),distance=int(d[i])))
 all_offsets['family_%d'%fi]=S.astype(np.uint8);candidate_total+=total
 row=dict(functionals=L.tolist(),offset_trials=total,non_torsion_pool=len(pool),survived_128_information_sets=len(S),strong_minimum_upper_estimate_histogram=dict(Counter(map(int,d))),history=history)
 # A remaining pool is useful only if four compatible nonzero cosets exist.
 if len(S)>=4:
  from itertools import combinations
  pairs=list(combinations(range(len(S)),2));DV=np.array([(S[i]-S[j])%8 for i,j in pairs]);pd,pc=restricted(DV,L,plans[160:]);adj=[0]*len(S)
  for (i,j),distance in zip(pairs,pd):
   if distance>=24:adj[i]|=1<<j;adj[j]|=1<<i
  clique=find_clique(adj,4)
  row['pair_differences_tested']=len(pairs);row['unverified_five_clique']=clique
  if clique:frontier.append(dict(functionals=L.tolist(),offsets=[np.zeros(24,dtype=int).tolist()]+S[clique].tolist(),status='unverified; exactpairdistancecertificate required'))
 rows.append(row);print('FAMILY',fi,'trials',total,'pool',len(pool),'strong survivors',len(S),'max upper',int(d.max()),flush=True)
receipt=dict(status='bounded seeded non-torsion maximin generation; all sampled minima are upper estimates',seed=SEED+5,offset_trials=candidate_total,rows=rows,frontier_candidates=frontier,strong_rejection_examples=badwitnesses,elapsed_seconds=time.monotonic()-start,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
np.savez_compressed(OUT/'index4_ascent_survivors.npz',**all_offsets);(OUT/'index4_ascent_results.json').write_text(json.dumps(receipt,indent=2)+'\n');print('SUMMARY',candidate_total,'frontier',len(frontier),'seconds',receipt['elapsed_seconds'],flush=True)
