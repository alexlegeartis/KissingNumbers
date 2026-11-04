"""Reconstruct both certificates without search or optimization."""
from pathlib import Path
import hashlib,json
import numpy as np
from verify import ROOT,sign_cube31
from verify_geometry import require

if __name__=='__main__':
    for record in json.loads((ROOT/'data/provenance.json').read_text())['records']:
        source=ROOT/record['baseline'];require(hashlib.sha256(source.read_bytes()).hexdigest()==record['baseline_sha256'],'Baseline hash')
        old=np.load(source,allow_pickle=False)['points'];delta=np.load(ROOT/record['delta'],allow_pickle=False)
        removed={row.tobytes() for row in delta['removed']};kept=np.array([row for row in old if row.tobytes() not in removed],np.int16)
        added=sign_cube31() if record['dimension']==31 else delta['added']
        points=np.unique(np.r_[kept,added],axis=0);expected=np.load(ROOT/record['certificate'],allow_pickle=False)['points']
        require(np.array_equal(points,expected),'Reconstructed exact coordinate array')
        print(json.dumps({'status':'PASS','dimension':record['dimension'],'points':len(points)}))
