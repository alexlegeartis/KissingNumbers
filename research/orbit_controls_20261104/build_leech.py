"""Build an exact Leech whole orbit from small algebraic generators."""
from pathlib import Path
import hashlib,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'common'))
from golay import golay24
from m24 import find_autos
from discovery.core import RunContext
from discovery.engines.exact_orbits import integer_orbit

def main():
 out=Path(__file__).resolve().parent;start=time.monotonic();G=golay24();words=set(map(tuple,G.tolist()))
 autos=find_autos(G,n_want=1,seed=20261104,tries=2)
 if not autos:raise ValueError('Bounded M24 generator search returned no automorphism')
 cycle=list(range(1,23))+[0,23]
 permutations=[cycle,autos[0].tolist()]
 generators=[]
 for p in permutations:
  if sorted(p)!=list(range(24)) or any(tuple(w[p]) not in words for w in G):raise ValueError('Golay permutation')
  M=np.eye(24,dtype=np.int64)[p]
  generators.append({'matrix':M.tolist(),'denominator':1})
 M=np.diag(1-2*G[1].astype(np.int64));generators.append({'matrix':M.tolist(),'denominator':1})
 # Derive the sextet containing a deterministic tetrad from all five Golay
 # octads containing it; every union of two sextet tetrads is an octad.
 tetrad={0,1,2,3};blocks=[tuple(sorted(tetrad))]
 for word in G[G.sum(axis=1)==8]:
  supp=set(map(int,np.flatnonzero(word)))
  if tetrad<=supp:blocks.append(tuple(sorted(supp-tetrad)))
 if len(blocks)!=6 or sorted(x for b in blocks for x in b)!=list(range(24)):raise ValueError('Golay sextet')
 for i in range(6):
  for j in range(i):
   indicator=[int(k in set(blocks[i]+blocks[j])) for k in range(24)]
   if tuple(indicator) not in words:raise ValueError('Sextet octad union')
 M=2*np.eye(24,dtype=np.int64)
 for block in blocks:
  for i in block:
   for j in block:M[i,j]=1-2*int(i==j)
 M[list(blocks[0])]*=-1
 generators.append({'matrix':M.tolist(),'denominator':2})
 seed=[4,4]+[0]*22
 program={'schema_version':1,'name':'Leech exact orbit control','nodes':[{'id':'sphere','op':'groups.integer_orbit','params':{'seed':seed,'generators':generators}}],'outputs':['sphere']}
 (out/'leech_program.json').write_text(json.dumps(program,indent=2)+'\n')
 context=RunContext(ROOT,{'max_materialized_points':250000,'max_orbit_points':250000,'max_rank_work':1000000000,'max_pairs':250000})
 print('generated four exact algebraic actions',flush=True)
 artifact=integer_orbit([],program['nodes'][0]['params'],context)
 print(json.dumps({'facts':artifact.facts,'proof':{k:artifact.verification.details[k] for k in ('seed_pairs_checked','unordered_pairs_covered','generator_images_checked')},'seconds':time.monotonic()-start}),flush=True)
 (out/'leech_artifact.json').write_text(json.dumps(artifact.to_dict(),separators=(',',':'))+'\n')
 (out/'leech_summary.json').write_text(json.dumps({'facts':artifact.facts,'verification':artifact.verification.to_dict(),'budget':context.budget,'program':'leech_program.json','artifact':'leech_artifact.json','seconds':time.monotonic()-start,'new_sota':False,'scope':'Known Leech control from generated Golay automorphism and exact sextet actions; not held-out discovery'},indent=2)+'\n')
if __name__=='__main__':main()
