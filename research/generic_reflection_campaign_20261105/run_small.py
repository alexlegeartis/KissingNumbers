"""Run one generic reflection operator on a supplied benchmark seed."""
import argparse,json,pathlib,sys,numpy as np
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]))
from discovery.core import RunContext
from discovery.engines.rational_reflection_search import search,quadratic_search,certify_weighted
from discovery.engines.projective_exchange import certify
from discovery.engines.quadratic_exchange import certify as qcertify
p=argparse.ArgumentParser();p.add_argument('dimension',type=int,choices=tuple(range(12,24)));p.add_argument('--axes',type=int,default=4);p.add_argument('--candidates',type=int,default=60);p.add_argument('--seconds',type=float,default=8);p.add_argument('--out',type=pathlib.Path);args=p.parse_args();where='data' if args.dimension in (19,20,21) else 'cache';z=np.load(f'research/generic_reflection_campaign_20261105/{where}/seed{args.dimension}.npz');rng=np.random.default_rng(11000+args.dimension);ctx=RunContext(pathlib.Path('.'),{'max_pairs':10000000000,'max_rank_work':250000000000,'max_materialized_points':500000,'max_codewords':300000,'max_gram_entries':3000000,'solver_seconds':30});common={'candidate_limit':args.candidates,'seconds':args.seconds}
if args.dimension in (13,18):
 A=z['A'].astype(int).tolist();B=z['B'].astype(int).tolist();idx=rng.choice(len(A),args.axes,replace=False);v=quadratic_search([],dict(common,source_a=A,source_b=B,axes_a=[A[int(i)] for i in idx],axes_b=[B[int(i)] for i in idx],metric=(z['g'].astype(int).tolist() if args.dimension==18 else [1]*13)),ctx);code=(certify_weighted if args.dimension==18 else qcertify)([v],{},ctx)
else:
 X=z['X'].astype(int).tolist();g=z['g'].astype(int).tolist();idx=rng.choice(len(X),args.axes,replace=False);v=search([],dict(common,source=X,axes=[X[int(i)] for i in idx],metric=g),ctx);code=certify([v],{},ctx)
detail=v.data['proposal'] if v.data.get('backend')=='weighted_quadratic_reflection' else v.data
result={'dimension':args.dimension,'source_count':len(A) if args.dimension in (13,18) else len(X),'verified_count':code.facts['cardinality'],'proposal_count':len(detail['selected_indices']),'axis_indices':idx.tolist(),'solver_status':detail['solver_status'],'verification':code.verification.status,'scope':'finite reflection pool; no automatic literature priority claim'};print(json.dumps(result,indent=2))
if args.out:
 if args.out.exists():raise FileExistsError(args.out)
 args.out.write_text(json.dumps({'summary':result,'proposal':v.to_dict(),'code':code.to_dict()}))
