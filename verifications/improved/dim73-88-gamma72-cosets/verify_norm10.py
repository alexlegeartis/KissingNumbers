"""Replay the explicit 73D addition and exact 77--79D existence proofs."""
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[2]/'common'))
from checkpoint import check_manifest,replay,temporary_source
from final import results

check_manifest(ROOT)
with temporary_source(ROOT) as source:
    replay(source,'dim73_surgery/verify_surgery.py')
    replay(source,'norm10/verify_shell10.py')
    replay(source,'norm10/norm10_bound_certificate.py')
    replay(source,'norm10/weighted_norm10_certificate.py')
    replay(source,'norm10/verify_weighted_independent.py')
    explicit=json.loads((source/'norm10/shell10_verification.json').read_text())
    whole=json.loads((source/'norm10/norm10_bound_certificate.json').read_text())
    fresh={73:explicit['bound']}|{int(d):r['bound'] for d,r in whole['dimensions'].items()}
    weighted=json.loads((source/'norm10/weighted_norm10_certificate.json').read_text())
    for dim,row in weighted['dimensions'].items():
        fresh[int(dim)]=max(fresh.get(int(dim),0),row['bound'])
    expected={r['dimension']:r['total'] for r in results() if r['dimension'] in fresh}
    assert fresh==expected
print(json.dumps({'status':'PASS','results':fresh}))
