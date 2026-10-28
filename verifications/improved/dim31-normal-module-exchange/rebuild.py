"""Rebuild the final point set from the pinned public PackingStar file and delta."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
ROOT=Path(__file__).resolve().parent


def main(baseline,output):
    manifest=json.loads((ROOT/'data/provenance.json').read_text(encoding='utf-8'))
    assert hashlib.sha256(baseline.read_bytes()).hexdigest()==manifest['packingstar_coordinate_sha256']
    with np.load(baseline,allow_pickle=False) as z:
        denominator=int(z['coordinate_denominator']);assert np.all(z['points']%denominator==0)
        x=(z['points']//denominator).astype(np.int8)
    with np.load(ROOT/'data/exchange_delta.npz',allow_pickle=False) as z:
        removed=z['removed'];added=z['added']
    remove={r.tobytes() for r in removed};keys={r.tobytes() for r in x}
    assert remove<=keys
    kept=x[[r.tobytes() not in remove for r in x]]
    assert not any(r.tobytes() in keys-remove for r in added)
    points=np.unique(np.r_[kept,added],axis=0)
    assert len(points)==len(x)-len(removed)+len(added)==251110
    digest=hashlib.sha256(points.tobytes()).hexdigest()
    assert digest==manifest['canonical_coordinate_sha256']
    np.savez_compressed(output,points=points,squared_norm=np.int32(32),inner_product_limit=np.int32(16),coordinate_denominator=np.int32(1))
    print(json.dumps(dict(points=len(points),removed=len(removed),added=len(added),canonical_coordinate_sha256=digest,output=str(output))))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('baseline',type=Path);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();main(a.baseline,a.output)
