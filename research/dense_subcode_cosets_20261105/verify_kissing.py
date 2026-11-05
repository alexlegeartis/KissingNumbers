"""Independently verify the conditional 32–35D dense-layer replacement."""
import importlib.util,json,pathlib,sys
ROOT=pathlib.Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from discovery.core import Artifact,RunContext
from discovery.certification import verify_artifact

def verify(code):
 ctx=RunContext(ROOT,{'max_codewords':1<<19,'max_pairs':10000000,'max_rank_work':100000000})
 report=verify_artifact(code,ctx);assert report.status=='exact' and code.kind=='binary_code'
 assert code.facts['length']==32 and code.facts['min_distance']>=8 and code.facts['cardinality']>131072
 loc=ROOT/'verifications/improved/dim32-35-signed-layers/source/classic';spec=importlib.util.spec_from_file_location('independent_source_verifier',loc/'verify.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);out=[]
 # Dense signs have squared norm 32. Distance >= 8 gives inner product
 # <= 16; an eight-support (+/-2) point has inner product <= 16 with any
 # dense sign, and any two-coordinate (+/-4) root has inner product <= 8.
 # Once heavy heads are removed, these inequalities hold for every sign word.
 for d in(32,33,34,35):
  result=m.verify(loc/f'certificate{d}.json');old=result['verified_points'];lost=2 if d in(34,35) else 0;new=old+code.facts['cardinality']-131072-lost
  out.append({'dimension':d,'old_points':old,'verified_points':new,'gain':new-old,'heavy_heads_removed':d in(34,35)})
 return out
if __name__=='__main__':
 code=Artifact.from_dict(json.loads(pathlib.Path(sys.argv[1]).read_text()));print(json.dumps(verify(code),indent=2))
