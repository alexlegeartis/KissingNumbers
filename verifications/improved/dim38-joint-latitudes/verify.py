"""Check complete integrity and replay both historical and current exact proofs."""
from pathlib import Path
import json
import shutil
import sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[2]/'common'))
from checkpoint import check_manifest,temporary_source,replay
from final import results

if __name__=='__main__':
    check_manifest(ROOT)
    expected=results()[0]
    with temporary_source(ROOT) as source:
        replay(source,'verify.py')
        old=json.loads((source/'verification.json').read_text())
        assert old['all_checks_passed'] and old['lower_bound']==610394
        refill=source.parent/'refill'
        shutil.copytree(ROOT/'refill',refill,ignore=shutil.ignore_patterns('__pycache__'))
        replay(refill,'verify.py')
        replay(refill,'verify_coordinates.py')
        replay(refill,'construct.py','--check')
        receipts=[json.loads((refill/name).read_text()) for name in ('verification.json','coordinate_verification.json')]
        assert all(r['all_checks_passed'] and r['lower_bound']==expected['total'] for r in receipts)
        assert receipts[0]['certificate_sha256']==receipts[1]['certificate_sha256']
    print(json.dumps(dict(status='PASS',results=[expected])))
