"""Reproducible coordinate-equality subcode campaign; each result is scoped."""
import argparse,json,pathlib,sys
import numpy as np
ROOT=pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from discovery.core import RunContext
from discovery.engines.binary_coset_discovery import search

p=argparse.ArgumentParser()
p.add_argument('--trials',type=int,default=1000)
p.add_argument('--start',type=int,default=0)
p.add_argument('--seconds',type=float,default=120)
p.add_argument('--max-candidates',type=int,default=6400)
p.add_argument('--seed',type=int,default=326913)
a=p.parse_args()
if a.start<0 or a.trials<a.start:p.error('require 0 <= start <= trials')
src=json.loads((ROOT/'verifications/improved/dim32-35-signed-layers/source/classic/certificate32.json').read_text())
rng=np.random.default_rng(a.seed)
for i in range(a.trials):
 block=sorted(map(int,rng.choice(32,4,replace=False)))
 if i<a.start:continue
 ctx=RunContext(ROOT,{'max_trials':1,'max_codewords':1<<18,'max_rank_work':32*(1<<18),'max_gram_entries':a.max_candidates**2,'solver_seconds':a.seconds})
 params={'length':32,'parent_generator':src['dense_generator'],'codimension':3,'minimum_distance':8,'trials':1,'seconds':a.seconds,'max_candidates':a.max_candidates,'coordinate_equalities':[block]}
 out=search([],params,ctx)
 print(json.dumps({'trial':i,'block':block,'search':out.data['trials'],'proposal':out.data['proposal']}),flush=True)
