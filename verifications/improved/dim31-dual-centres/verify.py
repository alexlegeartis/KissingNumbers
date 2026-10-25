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
        replay(source,'verify31.py','certificate31_241742.npz','--output','fresh.json')
        receipt=json.loads((source/'fresh.json').read_text())
        assert receipt['status']=='PASS EXACT' and receipt['points']==expected['total']
        replay(source,'materialize31.py','certificate31_241742.npz','fresh_points.npz')
        import numpy as np
        with np.load(source/'fresh_points.npz',allow_pickle=False) as fresh, np.load(source/'points31_241742.npz',allow_pickle=False) as recorded:
            assert set(fresh.files)==set(recorded.files)
            assert all(np.array_equal(fresh[k],recorded[k]) for k in fresh.files)
    print(json.dumps(dict(status='PASS',results=[expected])))
