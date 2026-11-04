"""Run the registered general exchange operator on the discovered 31D map.

The operator receives a source and rational isometry, with no output point
list, target count, cutoff or prescribed exchange. All cross conflicts are
constructed; source and output geometry are freshly independently replayed.
"""
from pathlib import Path
import gzip,hashlib,json,sys,time
from math import isqrt
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from discovery.core import Artifact,RunContext
from discovery.runner import registry

def main():
 out=Path(__file__).resolve().parent
 report=json.loads((out/'images31/report.json').read_text())
 source=ROOT/report['source'];raw=source.read_bytes()
 if hashlib.sha256(raw).hexdigest()!=report['source_sha256']:raise ValueError('Pinned source')
 action=next(r['action'] for r in report['records'] if r.get('gain',0)>0)
 if action['denominator']!=1:raise ValueError('Adapter expects this integer signed-permutation map')
 X=np.load(source,allow_pickle=False)['points']
 source_artifact=Artifact('spherical_code',{'coordinates':X.astype(np.int64).tolist(),'norm_squared':128})
 context=RunContext(ROOT,{'max_materialized_points':600000,'max_pairs':150000000000,
  'max_gram_entries':150000000,'max_rank_work':1000000000})
 started=time.monotonic();print('Registered selection.isometric_exchange: complete source replay, full cross graph, complete output replay',flush=True)
 result=registry()[0]['selection.isometric_exchange']([source_artifact],{'matrix':action['numerator']},context)
 with gzip.open(out/'engine31_artifact.json.gz','wt',encoding='utf8') as f:json.dump(result.to_dict(),f,separators=(',',':'));f.write('\n')
 norm=result.data['norm_squared'];scale=isqrt(128//norm)
 if norm*scale*scale!=128:raise ValueError('Exact portable coordinate scaling')
 points=np.array(result.data['coordinates'],np.int16)*scale
 np.savez_compressed(out/'engine31.npz',points=points,squared_norm=128,inner_product_limit=64,coordinate_denominator=2)
 proof=result.verification.details
 summary={'operator':'selection.isometric_exchange','facts':result.facts,
  'verification_status':result.verification.status,'seconds':time.monotonic()-started,
  'source':str(source.relative_to(ROOT)),'source_sha256':hashlib.sha256(raw).hexdigest(),
  'budget':context.budget,'isometric_exchange':{k:v for k,v in proof['isometric_exchange'].items() if k not in
   ('matching','left_indices','right_indices','cover_left','cover_right','right_source_indices','matrix_numerator')},
  'new_sota':False,'scope':'Seeded algebraic isometry search. Exact finite-union optimality and full configuration replay; no unrestricted optimum or worldwide-priority assertion.'}
 (out/'engine31_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':main()
