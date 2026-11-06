"""Replay inherited source, both new exact proofs and corruption controls."""
from pathlib import Path
import json,shutil,sys,tempfile
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[2]/'common'))
from checkpoint import check_manifest,replay
from final import results
if __name__=='__main__':
 check_manifest(ROOT);expected=results()[0]
 with tempfile.TemporaryDirectory(prefix='kissing38_hamming_') as d:
  d=Path(d)
  for name in ['source','refill','hamming']:shutil.copytree(ROOT/name,d/name,ignore=shutil.ignore_patterns('__pycache__'))
  replay(d/'source','verify.py');assert json.loads((d/'source/verification.json').read_text())['all_checks_passed']
  replay(d/'refill','verify.py');replay(d/'refill','verify_coordinates.py');replay(d/'refill','construct.py','--check')
  for name in ['verification.json','coordinate_verification.json']:
   r=json.loads((d/'refill'/name).read_text());assert r['all_checks_passed'] and r['lower_bound']==611522
  replay(d/'hamming','construct.py','--check');replay(d/'hamming','verify.py');replay(d/'hamming','verify_independent.py');replay(d/'hamming','negative_controls.py')
  a=json.loads((d/'hamming/verification.json').read_text());b=json.loads((d/'hamming/independent_verification.json').read_text());assert a['all_new_and_cross_layer_checks_passed'] and b['status']=='PASS' and a['lower_bound']==b['total']==expected['total'] and a['certificate_sha256']==b['certificate_sha256'];controls=json.loads((d/'hamming/negative_control_results.json').read_text());assert len(controls)==4 and all(r['rejected'] for r in controls)
 print(json.dumps(dict(status='PASS',results=[expected],proof_routes=2,negative_controls=4,inherited_input='Source Golay/Leech minimum-norm argument; no exhaustive Leech-pair enumeration.')))
