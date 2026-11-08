"""Run a relocated positive replay and four deliberately invalid certificates."""
from pathlib import Path
import tempfile,shutil,subprocess,sys,json,os
import numpy as np
HERE=Path(__file__).resolve().parent

def replay(path):
 p=subprocess.run([sys.executable,str(path/'verify.py')],cwd=path,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
 return p.returncode

def mutate(path,case):
 data=path/'data'
 if case in ('incumbent_hash','pool_hash'):
  p=data/('incumbent251200.npz' if case=='incumbent_hash' else 'repriced_pool.npz');z=np.load(p);d={k:z[k] for k in z.files};k='points' if case=='incumbent_hash' else 'candidates';d[k]=d[k].copy();d[k][0,0]+=1;np.savez_compressed(p,**d)
 elif case=='dual_sign':
  p=data/'lp_dual_certificate.npz';z=np.load(p);d={k:z[k] for k in z.files};d['dual_inequalities']=d['dual_inequalities'].copy();d['dual_inequalities'][0]=1;np.savez_compressed(p,**d)
 elif case=='clique_duplicate':
  p=data/'lp_constraints.json';d=json.loads(p.read_text());d['global_cliques'][0].append(d['global_cliques'][0][0]);p.write_text(json.dumps(d))
 else:raise ValueError(case)

records=[]
for case in ('positive','incumbent_hash','pool_hash','dual_sign','clique_duplicate'):
 with tempfile.TemporaryDirectory(prefix='joint31-exact-') as temp:
  p=Path(temp)/'package';shutil.copytree(HERE,p)
  if case!='positive':mutate(p,case)
  code=replay(p);assert (code==0)==(case=='positive'),(case,code)
  records.append({'case':case,'exit_code':code,'expected_acceptance':case=='positive'})
report={'status':'PASS','controls':records};(HERE/'negative_controls.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
