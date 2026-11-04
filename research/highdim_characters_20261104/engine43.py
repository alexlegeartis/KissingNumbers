"""Adapt the archived43D incumbent to the generic body/head search operator.

No proposed delta or target count is passed to the engine. The bodies, head
palette and complete-family envelopes retain source/imported ancestry.
"""
from pathlib import Path
import argparse,base64,gzip,hashlib,io,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
P=ROOT/'verifications/improved/dim41-46-p48-heads/partial43';sys.path.insert(0,str(P))
import verify_families as vf
from discovery.core import RunContext
from discovery.runner import run_program,output_artifact
OUT=Path(__file__).resolve().parent

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--stage2',action='store_true');args=ap.parse_args()
 start=time.monotonic();source=ROOT/'research/43d_head_exchange_20261104'
 raw=base64.b64decode((source/'baseline_2919782.npz.b64').read_bytes())
 if hashlib.sha256(raw).hexdigest()!='ff4042cfe8bc6bc79f0b33aa4259df41a08c013e37511e03b6833a18aff321c4':raise ValueError('archived incumbent hash')
 old=dict(np.load(io.BytesIO(raw),allow_pickle=False));g,res,mu,_=vf.geometry()
 palette=np.unique(old['head_numerators'],axis=0)
 if args.stage2:
  delta=json.loads((OUT/'engine43_delta.json').read_text());i=delta['remove']
  old['body_indices']=np.r_[np.delete(old['body_indices'],i),[x[0] for x in delta['add']]]
  old['head_numerators']=np.vstack((np.delete(old['head_numerators'],i,axis=0),delta['head_numerators']))
 X=np.load(P/'tiny_bodies.npz')['H'].astype(np.int64);q=np.load(P/'body_complement.npz')
 co=vf.mm(vf.mm(X,q['b'].T),vf.inverse(q['q'],3))
 if np.any(co%288):raise ValueError('coefficient residues')
 co//=288;lookup={row.tobytes():i for i,row in enumerate(res)}
 labels=[lookup[row.tobytes()] for row in co%3];model=dict(np.load(P/'configuration_43.npz'))
 beta=[]
 for i in range(6561):beta.append((10+model['r3']-18+mu[((g['digits'][model['classes']]-g['digits'][i])%3)@g['powers']]).tolist())
 params={'bodies':X.tolist(),'body_squared_scale':288,'body_norm_squared':10,
  'heads':palette.tolist(),'head_denominator':int(old['head_denominator']),
  'incumbent_bodies':old['body_indices'].tolist(),'incumbent_heads':old['head_numerators'].tolist(),
  'body_classes':labels,'fixed_double_body_bounds':beta,'fixed_body_norms':model['r3'].tolist(),
  'fixed_heads':model['head_numerators'].tolist(),'fixed_head_denominator':int(model['head_denominator']),
  'pure_heads':model['pure_directions'].tolist()}
 if args.stage2:params.update(max_removed=2,replacement_test_limit=500000)
 program={'schema_version':1,'name':'Finite body/head exchange on archived43D inputs',
  'description':'Seeded discovery replay through generic operator. No delta or target count supplied. Imported body/head/envelope data.',
  'nodes':[{'id':'exchange','op':'search.body_head_exchange','params':params}],'outputs':['exchange']}
 context=RunContext(ROOT,{'max_materialized_points':3000000,'max_codewords':12000000,
  'max_trials':500000 if args.stage2 else 512,'max_gram_entries':128000000,'max_rank_work':6000000000})
 print('running generic body/head product search',flush=True)
 receipt=run_program(program,context);artifact=output_artifact(receipt)
 stem='body_head43_stage2' if args.stage2 else 'body_head43'
 for name,data in [(stem+'_program.json.gz',program),(stem+'_receipt.json.gz',receipt)]:
  with gzip.open(OUT/name,'wt',encoding='utf8') as f:json.dump(data,f,separators=(',',':'));f.write('\n')
 proposal=artifact.data
 if artifact.facts['gain']>0:
  if len(proposal['added_body_indices'])!=len(proposal['removed_indices'])+1:raise ValueError('Exchange cardinality')
  if args.stage2:
   delta={k:proposal[k] for k in ('removed_indices','added_body_indices','added_head_numerators','head_denominator')}
  else:
   if len(proposal['removed_indices'])!=1:raise ValueError('Expected first replay shape changed; inspect before adapting')
   delta={'remove':proposal['removed_indices'][0],
    'add':list(map(list,zip(proposal['added_body_indices'],proposal['added_head_indices']))),
    'head_denominator':proposal['head_denominator'],'head_numerators':proposal['added_head_numerators']}
  (OUT/('engine43_stage2_delta.json' if args.stage2 else 'engine43_delta.json')).write_text(json.dumps(delta,indent=2)+'\n')
 summary={'facts':artifact.facts,'search':proposal['search'],'seconds':time.monotonic()-start,
  'verification_status':artifact.verification.status,'full_construction_replay':'Run iterate43.py, then replay.py --delta research/highdim_characters_20261104/engine43_combined_delta.json' if args.stage2 else 'Run replay.py --delta research/highdim_characters_20261104/engine43_delta.json',
  'scope':'Generic seeded discovery replay. Original discovery belongs to the separate43D task. No held-out credit, global optimality or worldwide-priority conclusion.'}
 (OUT/('engine43_stage2_summary.json' if args.stage2 else 'engine43_summary.json')).write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':main()
