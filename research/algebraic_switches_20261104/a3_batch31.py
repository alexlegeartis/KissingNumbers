"""Bounded adaptive A3-completion search near the newly filled root hole."""
from pathlib import Path
from itertools import combinations,product
import argparse,json,os,subprocess,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path);args=parser.parse_args()
 out=args.out or HERE/'a3_batch31';out.mkdir(exist_ok=False);source=HERE/'lifts31/engine_candidate.npz'
 X=np.load(source)['points'];known={row.tobytes() for row in X};cases=[]
 for pair in combinations([0,5,6,8,10,20],2):
  for sign,tail_sign in product((-1,1),repeat=2):
   r=np.zeros(31,np.int16);r[list(pair)]=[8,8*sign]
   s=np.zeros(31,np.int16);s[[24,25]]=[8,8*tail_sign]
   if all(v.tobytes() in known for v in (r,-r,s,-s)):cases.append((pair,sign,tail_sign))
 rng=np.random.default_rng(20261104);order=rng.permutation(len(cases)).tolist();records=[];start=time.monotonic()
 env=dict(os.environ);env['OPENBLAS_NUM_THREADS']='1'
 for index in order[:16]:
  if time.monotonic()-start>240:break
  pair,sign,tail_sign=cases[index];directory=out/f'trial{len(records):02}'
  command=[sys.executable,str(HERE/'sparse_lifts.py'),'31','--blockers','4','--core-pair',*map(str,pair),
   '--tail-pair','24','25','--core-pair-sign',str(sign),'--tail-pair-sign',str(tail_sign),
   '--source',str(source),'--out',str(directory),'--seconds','15']
  with (out/f'trial{len(records):02}.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,env=env)
  report=json.loads((directory/'report.json').read_text());gain=report['optimization']['gain']
  record={'directory':str(directory.relative_to(HERE)),'core_pair':list(pair),'core_sign':sign,
   'tail_pair':[24,25],'tail_sign':tail_sign,'candidates':report['candidate_pool'],
   'admitted':report['admitted_candidates'],'gain':gain,'seconds':report['seconds']}
  records.append(record);print(json.dumps(record),flush=True)
  (out/'report.json').write_text(json.dumps({'source':str(source.relative_to(ROOT)),'baseline':len(X),
   'available_cases':len(cases),'case_limit':16,'wall_budget':240,'records':records,
   'scope':'Generated A3 completions with retained orthogonal root pairs; finite point-level four-blocker selection. Independent full proof required for a gain.'},indent=2)+'\n')
  if gain>0:
   subprocess.run([sys.executable,str(HERE/'run_blocker_pool.py'),str(directory)],check=True,env=env)
   subprocess.run([sys.executable,str(ROOT/'verifications/improved/dim27-31-joint-refills/verify_geometry.py'),
    str(directory/'engine_candidate.npz'),'--report',str(directory/'full_geometry32.json')],check=True,env=env)
   break
if __name__=='__main__':main()
