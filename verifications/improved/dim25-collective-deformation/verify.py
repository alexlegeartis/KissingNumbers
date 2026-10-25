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
        replay(source,'verify.py','--out','fresh.json')
        receipt=json.loads((source/'fresh.json').read_text())
        assert receipt['status']=='EXACT VERIFICATION PASSED' and receipt['points']==expected['total']
    print(json.dumps(dict(status='PASS',results=[expected])))
