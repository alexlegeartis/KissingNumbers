"""Exact character kernels + joint nonlinear large-neighborhood line labels."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
from probe21 import D,line_ids,basis,comp,directions,seed
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from discovery.core import RunContext
from discovery.engines.characters import family_problem,optimize_families,check_assignment
OUT=Path(__file__).resolve().parent

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--rounds',type=int,default=12);ap.add_argument('--seconds',type=float,default=4);args=ap.parse_args()
 g=[[(D[m]>>j)&1 for j in range(21)] for m in basis]
 shapes=[[2 if D[m]>>j&1 else 3 for j in range(21)] for m in line_ids]
 q=np.eye(10,dtype=int).tolist();context=RunContext(ROOT,{'max_rank_work':100000000,'max_gram_entries':3000000})
 p=family_problem(g,shapes,None,q,context);table=p['forbidden'];expected=np.zeros_like(table)
 for label,d in enumerate(directions):
  for a in range(21):
   for b in range(21):expected[d,a,b]=(a==label or b==label or a==b)
 expected[0]=True
 if not np.array_equal(table,expected):raise ValueError('Walsh kernel differs from independently derived line rule')
 assignment=seed[:];rng=np.random.default_rng(20261104);history=[];started=time.monotonic()
 for r in range(args.rounds):
  holes=np.flatnonzero(np.array(assignment)<0);v=int(rng.choice(holes))
  # Includes all first and second Cayley neighbors of the candidate hole.
  active={v}
  for d in directions:
   active.add(v^d)
   for e in directions:active.add(v^d^e)
  if r%3==2:
   # Break a larger coupled block rather than only one hole neighborhood.
   w=int(rng.choice(holes));active.update(w^d for d in directions);active.add(w)
  assignment,diagnostic=optimize_families(p,context,args.seconds,assignment,sorted(active))
  row={'round':r,'center':v,'selected':check_assignment(p,assignment),'diagnostic':diagnostic,'seconds':time.monotonic()-started}
  print(json.dumps(row),flush=True);history.append(row)
 out={'history':history,'best_assignment':assignment,'component_size':1024,'selected_per_component':check_assignment(p,assignment),'four_component_added_count':4*check_assignment(p,assignment),'verification_status':'exact finite compatibility only; solver diagnostics numerical','new_sota':False,'scope':'Nonlinear whole component line labels; fixed odd-octad/root base. Finite time-limited neighborhoods do not bound all labels.'}
 (OUT/'discrete21.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
