"""Joint base/cap orbit assembly with new mixed-amplitude norm72 templates.

Avoids the parent's flat-coset/norm-six-tier catalogue: new shapes must have
more than eight nonzero head entries and at least two head magnitudes.
"""
from pathlib import Path
from collections import Counter
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from discovery.core import RunContext
from discovery.engines.characters import assemble_orbits

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--seed',type=int,default=20261104);ap.add_argument('--candidates',type=int,default=192);ap.add_argument('--seconds',type=float,default=30);args=ap.parse_args()
 rng=np.random.default_rng(args.seed);out=Path(__file__).resolve().parent;start=time.monotonic()
 old=json.load(open(ROOT/'research/orbit_controls_20261104/control17_artifact.json'))['data']['coordinates']
 g=[[1]*16+[0]]+[[x>>j&1 for x in range(16)]+[0] for j in range(4)]
 signs=np.array([[1-2*(sum((m>>j&1)*g[j][x] for j in range(5))%2) for x in range(17)] for m in range(32)],dtype=np.int64)
 metric=np.array([1]*16+[2],np.int64)
 def key(row):return min(tuple(map(int,r)) for r in np.asarray(row)*signs)
 groups={}
 for row in old:groups.setdefault(key(row),set()).add(tuple(row))
 originals=sorted(groups)
 for row in originals:
  if set(map(tuple,signs*np.array(row)))!=groups[row]:raise ValueError('Incumbent not invariant under RM1')
 print(json.dumps({'incumbent_orbits':len(originals),'incumbent_points':len(old),'orbit_sizes':dict(Counter(map(len,groups.values())))}),flush=True)
 # Uniform bounded integer-square compositions; no record count or coordinates
 # determine the new candidate amplitudes.
 dp=[[0]*73 for _ in range(17)];dp[0][0]=1
 for k in range(1,17):
  for s in range(73):dp[k][s]=sum(dp[k-1][s-a*a] for a in range(5) if a*a<=s)
 def sample(total):
  vals=[]
  for k in range(16,0,-1):
   choices=[a for a in range(5) if a*a<=total and dp[k-1][total-a*a]]
   probabilities=np.array([dp[k-1][total-a*a] for a in choices],float);probabilities/=probabilities.sum()
   a=int(rng.choice(choices,p=probabilities));vals.append(a);total-=a*a
  return np.array(vals)*rng.choice([-1,1],16)
 census=Counter();proposals={};attempts=0;target=max(1200,args.candidates*6)
 while len(proposals)<target and attempts<100000:
  attempts+=1;h=int(rng.choice([-3,-2,-1,0,1,2,3]));head=sample(72-2*h*h)
  if np.count_nonzero(head)<=8 or len(set(abs(head[head!=0])))<2:continue
  row=np.r_[head,h];dots=(signs*row)@(metric*row)
  if np.any((2*dots>72)&(dots!=72)):continue
  k=key(row)
  if k in groups or k in proposals:continue
  # Exact blocker pricing is diagnostic, not a fixed-base admission test.
  dots=(signs*row)@(np.array(originals)*metric).T
  blockers=np.flatnonzero(dots.max(axis=0)>36)
  cost=sum(len(groups[originals[int(i)]]) for i in blockers)
  proposals[k]=(cost,h)
  census[(h,cost)]+=1
 print(json.dumps({'generated':len(proposals),'attempts':attempts,'minimum_blocker_cost':min(v[0] for v in proposals.values())}),flush=True)
 cheapest=sorted(proposals,key=lambda k:(proposals[k][0],k))[:args.candidates//2]
 rest=sorted(set(proposals)-set(cheapest));rng.shuffle(rest)
 selected=cheapest+rest[:args.candidates-len(cheapest)]
 shapes=[list(x) for x in originals+selected]
 context=RunContext(ROOT,{'max_gram_entries':len(shapes)**2*32+1,'max_rank_work':len(shapes)**2*200+1,'max_pairs':1000000000,'max_materialized_points':100000,'solver_seconds':max(30,args.seconds)})
 params={'generators':g,'shapes':shapes,'weights':metric.tolist(),'initial_shapes':list(range(len(originals))),'seconds':args.seconds}
 program={'schema_version':1,'name':'17D mixed-amplitude coupled orbit search','nodes':[{'id':'sphere','op':'search.character_orbits','params':params}],'outputs':['sphere']}
 artifact=assemble_orbits([],params,context)
 stem=f'coupled17-s{args.seed}'
 (out/f'{stem}-program.json').write_text(json.dumps(program,separators=(',',':'))+'\n')
 (out/f'{stem}-artifact.json').write_text(json.dumps(artifact.to_dict(),separators=(',',':'))+'\n')
 report={'seed':args.seed,'candidate_count':len(selected),'generated_candidates':len(proposals),'attempts':attempts,'census':[{'height':h,'old_blocker_points':cost,'candidate_orbits':n} for (h,cost),n in sorted(census.items())], 'facts':artifact.facts,'search':artifact.data['search'],'seconds':time.monotonic()-start,'budget':context.budget,'new_sota':False,'scope':'Known5730 control plus sampled mixed-amplitude RM1 orbits. Every base/cap orbit removable. No full alphabet or real-coordinate exclusion.'}
 (out/f'{stem}-summary.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({'facts':artifact.facts,'search':artifact.data['search'],'seconds':report['seconds']}),flush=True)
if __name__=='__main__':main()
