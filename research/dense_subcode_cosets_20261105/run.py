"""Discover a larger dense binary layer by coset recombination."""
import argparse,json,pathlib,sys,importlib.util
ROOT=pathlib.Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from discovery.core import RunContext
from discovery.engines.binary_coset_discovery import search,certify
p=argparse.ArgumentParser();p.add_argument('--codimension',type=int,choices=(3,4),default=3);p.add_argument('--trials',type=int,default=10);p.add_argument('--seconds',type=float,default=60);p.add_argument('--max-candidates',type=int,default=4000);p.add_argument('--output',type=pathlib.Path);args=p.parse_args()
loc=ROOT/'verifications/improved/dim32-35-signed-layers/source/classic';src=json.loads((loc/'certificate32.json').read_text());states=1<<(15+args.codimension);ctx=RunContext(ROOT,{'max_trials':max(args.trials,1),'max_codewords':states,'max_rank_work':max(50000000,args.trials*32*states),'max_gram_entries':args.max_candidates**2,'solver_seconds':max(args.seconds,1)})
r=search([] ,{'length':32,'parent_generator':src['dense_generator'],'codimension':args.codimension,'minimum_distance':8,'trials':args.trials,'seconds':args.seconds,'max_candidates':args.max_candidates,'seed':323333,'scramble_steps':80},ctx)
print(json.dumps({'trials':r.data['trials'],'proposal_found':r.data['proposal'] is not None},indent=2))
if r.data['proposal']:
 code=certify([r],{},ctx);print(json.dumps({'exact_binary_code':code.facts},indent=2))
 import importlib.util
 v=importlib.util.spec_from_file_location('dense_kissing_verifier',pathlib.Path(__file__).with_name('verify_kissing.py'));m=importlib.util.module_from_spec(v);v.loader.exec_module(m);print(json.dumps({'verified_kissing_counts':m.verify(code)},indent=2))
 if args.output:args.output.write_text(json.dumps(code.to_dict()))
