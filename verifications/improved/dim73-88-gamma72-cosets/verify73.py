"""Replay the 73D non-lattice surgery in an isolated temporary copy.

Use --rebuild to regenerate the Leech lists and targeted completion as well.
"""
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[2]/'common'))
from checkpoint import check_manifest,replay,temporary_source
from final import results

if __name__=='__main__':
    check_manifest(ROOT)
    with temporary_source(ROOT) as source:
        if '--rebuild' in sys.argv:
            replay(source,'dim73_surgery/lattice/verify.py','--rebuild')
        replay(source,'dim73_surgery/verify_surgery.py')
        report=json.loads((source/'dim73_surgery/surgery_verification.json').read_text())
        assert next(r['total'] for r in results() if r['dimension']==73)==report['bound']
    print(json.dumps(dict(status='PASS',dimension=73,bound=report['bound'])))
