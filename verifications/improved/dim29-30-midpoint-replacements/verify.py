"""Check integrity and replay both independent midpoint proofs in isolation."""
from pathlib import Path
import json
import sys
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[2]/'common'))
from checkpoint import check_manifest,temporary_source,replay
from final import results

if __name__ == '__main__':
    check_manifest(ROOT)
    expected=results()
    with temporary_source(ROOT) as source:
        for row in expected:
            name='certificate%d_%d.npz'%(row['dimension'],row['total'])
            output='fresh%d.json'%row['dimension']
            replay(source,'verify.py',name,'--output',output)
            receipt=json.loads((source/output).read_text())
            assert receipt['status']=='PASS EXACT' and receipt['points']==row['total']
            expanded='fresh_points%d.npz'%row['dimension']
            replay(source,'materialize.py',name,expanded)
            with np.load(source/expanded,allow_pickle=False) as fresh, np.load(
                    source/('points%d_%d.npz'%(row['dimension'],row['total'])),
                    allow_pickle=False) as recorded:
                assert set(fresh.files)==set(recorded.files)=={'a','b'}
                assert all(np.array_equal(fresh[key],recorded[key]) for key in fresh.files)
    print(json.dumps(dict(status='PASS',results=expected)))
