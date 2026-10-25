"""Check complete package integrity and replay the proof in a temporary directory."""
from pathlib import Path
import json, shutil, sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[2]/'common'))
from checkpoint import check_manifest,temporary_source,replay
from final import results

if __name__=='__main__':
    check_manifest(ROOT)
    expected=results()[0]
    with temporary_source(ROOT) as source:
        replay(source,'verify.py')
        receipt=json.loads((source/'verification.json').read_text())
        assert receipt['all_checks_passed'] and receipt['lower_bound']==expected['total']
    print(json.dumps(dict(status='PASS',results=[expected])))
