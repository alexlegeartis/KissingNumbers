"""Discover polynomial-image supports, then search partial signed-line exchanges."""
import argparse,json,sys,time
from pathlib import Path
p=Path(__file__).resolve().parent;root=p.parents[1];sys.path.insert(0,str(root))
from discovery.core import RunContext
from discovery.engines.finite_field_search import generate,exchange
parser=argparse.ArgumentParser();parser.add_argument('--support-limit',type=int,default=11);parser.add_argument('--seconds',type=int,default=120);parser.add_argument('--output',default=str(p/'proposal.json'));a=parser.parse_args();c=json.load(open(p/'seed32.json'));ctx=RunContext(root,{'max_trials':120,'max_rank_work':150000000,'solver_seconds':max(1,a.seconds)})
params={'m':5,'modulus':37,'flat_rank':3,'exponents':[3,5,1],'coefficient_options':[[1],[1,2,3],list(range(32))], 'seed_supports':c['full_even_supports'],'max_seed_conflicts':3}
t=time.time();pool=generate([],params,ctx);print('polynomials',pool.data['polynomials_examined'],'permutations',pool.data['permutations'],'unique images',pool.data['unique_images'],'candidate supports',len(pool.data['supports']),'seconds',round(time.time()-t,1),flush=True)
if not pool.data['supports']:raise SystemExit('No candidate supports')
r=exchange([pool],{'seed_supports':c['full_even_supports'],'exceptional_lines':c['exceptional_antipodal_lines'],'parities':[0,1],'candidate_limit':min(a.support_limit,len(pool.data['supports'])),'seconds':a.seconds},ctx)
json.dump(r.to_dict(),open(a.output,'w'));print('proposal',r.facts,'solver',r.data['solver_status'],'output',a.output,flush=True)
