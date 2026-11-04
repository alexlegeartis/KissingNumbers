"""Bounded repeated engine search; final combined delta needs full exact replay."""
from pathlib import Path
import argparse,base64,gzip,io,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from discovery.core import RunContext
from discovery.runner import run_program,output_artifact
OUT=Path(__file__).resolve().parent

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--extra-rounds',type=int,default=3);args=ap.parse_args()
 with gzip.open(OUT/'body_head43_program.json.gz','rt') as f:program=json.load(f)
 params=program['nodes'][0]['params'];original=[(i,int(b),h) for i,(b,h) in enumerate(zip(params['incumbent_bodies'],params['incumbent_heads']))];records=original[:]
 first=json.loads((OUT/'engine43_delta.json').read_text());second=json.loads((OUT/'engine43_stage2_delta.json').read_text())
 def apply(removed,bodies,heads):
  nonlocal records
  removed=set(removed);records=[r for i,r in enumerate(records) if i not in removed]+[(None,b,h) for b,h in zip(bodies,heads)]
 apply([first['remove']],[x[0] for x in first['add']],first['head_numerators'])
 apply(second['removed_indices'],second['added_body_indices'],second['added_head_numerators'])
 params.update(max_removed=2,replacement_test_limit=500000)
 context=RunContext(ROOT,{'max_materialized_points':3000000,'max_codewords':12000000,'max_trials':500000,'max_gram_entries':128000000,'max_rank_work':6000000000})
 history=[];start=time.monotonic()
 for round in range(3,3+args.extra_rounds):
  params['incumbent_bodies']=[b for _,b,h in records];params['incumbent_heads']=[h for _,b,h in records]
  receipt=run_program(program,context);artifact=output_artifact(receipt);d=artifact.data
  with gzip.open(OUT/f'body_head43_round{round}_receipt.json.gz','wt') as f:json.dump(receipt,f,separators=(',',':'));f.write('\n')
  report={'round':round,'facts':artifact.facts,'search':d['search'],'seconds':time.monotonic()-start};history.append(report);print(json.dumps(report),flush=True)
  if artifact.facts['gain']<=0:break
  apply(d['removed_indices'],d['added_body_indices'],d['added_head_numerators'])
 kept={i for i,b,h in records if i is not None};added=[(b,h) for i,b,h in records if i is None]
 combined={'removed_indices':[i for i,b,h in original if i not in kept],'added_body_indices':[b for b,h in added],
  'head_numerators':[h for b,h in added],'head_denominator':params['head_denominator']}
 if len(records)!=len(original)+len(added)-len(combined['removed_indices']):raise ValueError('Combined delta count')
 (OUT/'engine43_combined_delta.json').write_text(json.dumps(combined,indent=2)+'\n')
 (OUT/'engine43_iterations.json').write_text(json.dumps({'history':history,'initial_tiny_points':len(original),'final_tiny_points':len(records),'gain':len(records)-len(original),'full_replay_required':True,'new_sota':False},indent=2)+'\n')
 print('combined tiny points',len(records),'gain',len(records)-len(original),flush=True)
if __name__=='__main__':main()
