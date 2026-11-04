"""Run a portable registered exchange DAG on an exact geometric point pool."""
from pathlib import Path
import argparse,hashlib,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from discovery.core import RunContext
from discovery.runner import run_program,output_artifact

def main():
 ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path);args=ap.parse_args();p=args.directory
 report=json.loads((p/'report.json').read_text());source=Path(report['source'])
 if not source.is_absolute():source=ROOT/source
 if hashlib.sha256(source.read_bytes()).hexdigest()!=report['source_sha256']:raise ValueError('Source hash')
 X=np.load(source)['points'].astype(np.int64);pool=np.load(p/'pool.npz')
 Z=pool['candidates'].astype(np.int64);ptr=pool['blocker_indptr'];edges=pool['blocker_indices']
 if not np.all(np.sum(Z*Z,axis=1)==128):raise ValueError('Candidate norms')
 blockers=[edges[ptr[i]:ptr[i+1]].tolist() for i in range(len(Z))]
 conflict=np.argwhere(np.triu(Z@Z.T>64,1)).tolist()
 params={'blockers':blockers,'candidate_conflicts':conflict,'seconds':120}
 program={'schema_version':1,'name':'Partial signed-orbit shared-deletion selection',
  'provenance':{'pool_sha256':hashlib.sha256((p/'pool.npz').read_bytes()).hexdigest()},
  'nodes':[{'id':'exchange','op':'selection.blocker_exchange','params':params}],'outputs':['exchange']}
 context=RunContext(ROOT,{'solver_seconds':120})
 receipt=run_program(program,context);proposal=output_artifact(receipt)
 (p/'engine_program.json').write_text(json.dumps(program,indent=2)+'\n')
 (p/'engine_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 selected=proposal.data['selected_candidate_indices'];removed=proposal.data['removed_source_indices']
 keep=np.ones(len(X),bool);keep[removed]=False;final=np.unique(np.r_[X[keep],Z[selected]],axis=0)
 if len(final)!=len(X)+proposal.facts['gain']:raise ValueError('Actual cardinality')
 lookup={row.tobytes():i for i,row in enumerate(final)}
 maximum=-128
 for row in Z[selected]:
  values=np.zeros(len(final),np.int64)
  for j in np.flatnonzero(row):values+=final[:,j]*int(row[j])
  values[lookup[row.tobytes()]]=-128
  if values.max()>64:raise ValueError('Direct integer changed-pair condition')
  maximum=max(maximum,int(values.max()))
 np.savez_compressed(p/'engine_candidate.npz',points=final.astype(np.int16),squared_norm=128,
                     inner_product_limit=64,coordinate_denominator=2)
 summary={'operator':'selection.blocker_exchange','baseline':len(X),'points':len(final),
  'facts':proposal.facts,'selected':selected,'removed':removed,'direct_integer_maximum':maximum,
  'changed_pairs':len(selected)*(len(final)-1),'verification_status':proposal.verification.status,
  'scope':'Registered finite-oracle selection plus direct exact changed-pair checks; full geometry replay still required'}
 (p/'engine_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
