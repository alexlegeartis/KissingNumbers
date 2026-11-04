"""Seeded five-quarter-coset candidate generation with rank-two parity kernels.
Every sampled minimum is an UPPER estimate, never a lower-distance certificate.
"""
from pathlib import Path
from itertools import combinations
from collections import Counter
import json,time,hashlib
import numpy as np
from search import info_set,canonical,H,SEED,RAW
OUT=Path(__file__).resolve().parent


def functional_key(L):
 A=L.copy()%2;at=0
 for j in range(12):
  ids=np.flatnonzero(A[at:,j])
  if not len(ids):continue
  i=at+int(ids[0]);A[[at,i]]=A[[i,at]]
  for k in range(2):
   if k!=at and A[k,j]:A[k]^=A[at]
  at+=1
  if at==2:return tuple(map(tuple,A.tolist()))
 return None


def plans_for(L,infos,rng):
 result=[]
 for I,A,R in infos:
  U=L@R.T%2;encoded=U[0]+2*U[1]
  pairs=[(j,k)for j in range(12)for k in range(j)if encoded[j] and encoded[k]and encoded[j]!=encoded[k]]
  j,k=pairs[int(rng.integers(len(pairs)))];solve=np.zeros((4,2),dtype=np.int64)
  for a in range(2):
   for b in range(2):solve[int(a*encoded[j])^int(b*encoded[k])]=[a,b]
  result.append((I,A,R,j,k,solve))
 return result


def restricted(V,L,plans):
 best=np.full(len(V),97,dtype=np.int16);messages=np.zeros((len(V),12),dtype=np.int16)
 for I,A,R,j,k,solve in plans:
  chosen=V[:,I];coeff=chosen@R%8;syndrome=coeff@L.T%2;label=syndrome[:,0]+2*syndrome[:,1]
  bits=solve[label];residual=(V-chosen@A)%8
  for evenj in range(0,8,2):
   dj=bits[:,0]+evenj
   for evenk in range(0,8,2):
    dk=bits[:,1]+evenk
    W=(residual-dj[:,None]*A[j]-dk[:,None]*A[k])%8;d=H[W].sum(axis=1)
    improve=d<best
    if np.any(improve):
     best[improve]=d[improve]
     messages[improve]=(coeff[improve]+dj[improve,None]*R[j]+dk[improve,None]*R[k])%8
 if np.any(messages@L.T%2):raise ValueError('oracle returned message outside rank-two kernel')
 return best,messages


def quotient_dedup(V,L,base):
 I,A,R=base;syndrome=canonical(V,base);parity=(V[:,I]@R@L.T)%2
 keys=np.column_stack([syndrome,parity]);_,indices=np.unique(keys,axis=0,return_index=True)
 return V[np.sort(indices)]


def find_clique(adjacency,target):
 def visit(candidates,chosen):
  if len(chosen)==target:return chosen
  while candidates:
   if candidates.bit_count()<target-len(chosen):return None
   bit=candidates&-candidates;v=bit.bit_length()-1;candidates^=bit
   found=visit(candidates&adjacency[v],chosen+[v])
   if found is not None:return found
  return None
 return visit((1<<len(adjacency))-1,[])


def main():
 start=time.monotonic();rng=np.random.default_rng(SEED+4);G=np.array(json.loads(RAW.read_text())['generator'],dtype=np.int64)
 base=info_set(G,range(24));outside=[j for j in range(24)if j not in base[0]]
 infos=[info_set(G,rng.permutation(24))for _ in range(640)]
 families=[np.vstack([G[:,i]%2,G[:,(i+7)%24]%2])for i in range(16)]
 families += [rng.integers(0,2,(2,12))for _ in range(16)]
 seen=set();rows=[];frontier=[];bad_witnesses=[];offset_total=0;pair_total=0;clique_total=0;all_vertex_data={}
 for family_index,L in enumerate(families):
  key=functional_key(L)
  if key is None or key in seen:continue
  seen.add(key);L=np.array(key,dtype=np.int64);plans=plans_for(L,infos,rng)
  encoded=L[0]+2*L[1];p,q=next((p,q)for p in range(12)for q in range(p)if encoded[p]and encoded[q]and encoded[p]!=encoded[q])
  E=np.array([np.zeros(24,dtype=np.int64),G[p],G[q],(G[p]+G[q])%8])
  random=np.zeros((256,24),dtype=np.int64);random[:,outside]=rng.integers(0,8,(256,12))
  even=np.zeros((128,24),dtype=np.int64);even[:,outside]=2*rng.integers(0,4,(128,12))
  torsion=np.zeros((128,24),dtype=np.int64);torsion[:,outside]=4*rng.integers(0,2,(128,12))
  sparse=[4*np.eye(24,dtype=np.int64)[i]for i in range(24)]
  sparse += [4*(np.eye(24,dtype=np.int64)[i]+np.eye(24,dtype=np.int64)[j])for i,j in combinations(range(24),2)]
  bits=((np.arange(24)[:,None]>>np.arange(5))&1).T;F=np.vstack([np.ones(24,dtype=int),bits,*[bits[i]*bits[j]for i in range(5)for j in range(i)]])
  quadratic=canonical(rng.integers(0,8,(128,len(F)))@F%8,base)
  offsets=canonical(np.vstack([random,even,torsion,sparse,quadratic,np.zeros((1,24),dtype=int)]),base)
  V=quotient_dedup(np.vstack([(offsets+e)%8 for e in E]),L,base);V=V[np.any(V,axis=1)]
  d,c=restricted(V,L,plans[:8]);offset_total+=len(V);eligible=np.flatnonzero(d>=24)
  # Bounded graph: keep controls, high upper estimates, and a seeded diverse sample.
  controls=[int(np.flatnonzero(np.all(V==e,axis=1))[0])for e in E[1:]]
  preferred=list(np.argsort(d)[-64:]);rest=[i for i in eligible if i not in preferred and i not in controls]
  selected=sorted(set(controls+preferred+list(rng.choice(rest,min(32,len(rest)),replace=False))))
  selected=[i for i in selected if d[i]>=24];S=V[selected];pairs=list(combinations(range(len(S)),2));pair_total+=len(pairs)
  differences=np.array([(S[i]-S[j])%8 for i,j in pairs],dtype=np.int64)
  pd,pc=restricted(differences,L,plans[8:16]) if len(pairs)else(np.array([],dtype=int),np.zeros((0,12),dtype=int))
  adjacency=[0]*len(S)
  for (i,j),distance in zip(pairs,pd):
   if distance>=24:adjacency[i]|=1<<j;adjacency[j]|=1<<i
  vertices= len(S);weak_edges=sum(a.bit_count()for a in adjacency)//2;attempts=0;edge_cuts=0;vertex_cuts=0;status='no_five_clique_in_retained_upper_graph';best_upper_clique=4
  # Adding D itself reduces the target to a four-clique among these vertices.
  for attempt in range(96):
   clique=find_clique(adjacency,4)
   if clique is None:break
   attempts+=1;clique_total+=1;best_upper_clique=5
   offsets5=np.vstack([np.zeros((1,24),dtype=np.int64),S[clique]])
   ten=list(combinations(range(5),2));DV=np.array([(offsets5[i]-offsets5[j])%8 for i,j in ten])
   strong,sm=restricted(DV,L,plans[16:144]);bad=np.flatnonzero(strong<24)
   if not len(bad):strong,sm=restricted(DV,L,plans[144:]);bad=np.flatnonzero(strong<24)
   if not len(bad):
    frontier.append(dict(functionals=L.tolist(),offsets=offsets5.tolist(),sampled_minimum_upper_estimates=strong.tolist(),information_set_queries=624,status='unverified five-coset candidate; exact distance proof required'))
    status='unverified_frontier';break
   for index in bad:
    i,j=ten[int(index)];distance=int(strong[index]);message=sm[index].tolist()
    if distance<=0:raise ValueError('quotient dedup failed: coincident D cosets')
    bad_witnesses.append(dict(functionals=L.tolist(),offset_difference=DV[index].tolist(),parent_message=message,distance=distance))
    if i==0:
     v=clique[j-1]
     for a in range(vertices):adjacency[a]&=~(1<<v)
     adjacency[v]=0;vertex_cuts+=1
    else:
     a,b=clique[i-1],clique[j-1];adjacency[a]&=~(1<<b);adjacency[b]&=~(1<<a);edge_cuts+=1
  else:status='clique_refinement_budget_cutoff'
  row=dict(family_index=family_index,functionals=L.tolist(),offsets_tested=len(V),survived_first_upper_estimate=int(len(eligible)),retained_vertices=vertices,weak_upper_graph_edges=weak_edges,five_cliques_refined=attempts,strong_edge_cuts=edge_cuts,strong_vertex_cuts=vertex_cuts,status=status,best_initial_minimum_upper_estimate=int(d.max()),best_found_upper_graph_clique_size=best_upper_clique)
  rows.append(row);all_vertex_data['family_%d_vertices'%family_index]=S.astype(np.uint8)
  print(row,flush=True)
 receipt=dict(status='bounded candidate generation; sampled minima are upper estimates only',seed=SEED+4,
  forecast=dict(subcode_index=4,quarter_coset_cardinality=2**34,required_even_cosets=5,odd_parent_cardinality=2**36,total_if_certified=9*2**34,repository_comparison=142162362560),
  functionals_tested=len(rows),offsets_tested=offset_total,pair_differences_tested=pair_total,five_cliques_refined=clique_total,frontier_candidates=frontier,rows=rows,strong_rejection_witnesses=bad_witnesses,elapsed_seconds=time.monotonic()-start,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 np.savez_compressed(OUT/'index4_vertices.npz',**all_vertex_data)
 (OUT/'index4_results.json').write_text(json.dumps(receipt,indent=2)+'\n');print('SUMMARY',len(rows),offset_total,pair_total,clique_total,'frontier',len(frontier),'seconds',receipt['elapsed_seconds'],flush=True)
if __name__=='__main__':main()
