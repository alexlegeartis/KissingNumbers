"""Finite Q(sqrt(3)) exchange search from the 1154-point 13D code."""
import argparse,json,pathlib,sys
import numpy as np
ROOT=pathlib.Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from discovery.core import RunContext
from discovery.engines.quadratic_exchange import search,certify
p=argparse.ArgumentParser();p.add_argument('--seconds',type=float,default=120);p.add_argument('--output',type=pathlib.Path);a=p.parse_args()
z=np.load(pathlib.Path(__file__).with_name('pool.npz'));X=z['source'].astype(int);Q=z['candidates'].astype(int)
assert X.shape==(1154,26) and Q.shape==(1728,26)
ctx=RunContext(ROOT,{'max_materialized_points':10000,'max_pairs':10000000,'max_rank_work':100000000,'solver_seconds':max(a.seconds,1)})
r=search([],{'source_a':X[:,:13].tolist(),'source_b':X[:,13:].tolist(),'candidate_a':Q[:,:13].tolist(),'candidate_b':Q[:,13:].tolist(),'minimum_gain':1,'seconds':a.seconds},ctx)
c=certify([r],{},ctx)
print(json.dumps({'source_points':len(X),'candidates':len(Q),'solver_status':r.data['solver_status'],'selected':len(r.data['selected_indices']),'removed':len(r.data['removed_indices']),'exact_facts':c.facts},indent=2))
if a.output and r.data['selected_indices']:a.output.write_text(json.dumps(r.data))
