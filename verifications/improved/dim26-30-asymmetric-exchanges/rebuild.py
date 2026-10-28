"""Apply the exact exchange to the pinned public input; verify the final set."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent


def main(d,source,output):
    meta=json.loads((ROOT/'data/provenance.json').read_text());record=next(r for r in meta['records'] if r['dimension']==d)
    assert hashlib.sha256(source.read_bytes()).hexdigest()==record['source_sha256']
    with np.load(source,allow_pickle=False) as z:
        den=int(z['coordinate_denominator']);assert np.all(z['points']%den==0);x=(z['points']//den).astype('<i2')
    with np.load(ROOT/f'data/exchange{d}.npz',allow_pickle=False) as z:removed=z['removed'];added=z['added']
    keys=lambda a:np.ascontiguousarray(a).view(np.dtype((np.void,2*d))).ravel()
    delete=np.isin(keys(x),keys(removed));assert int(delete.sum())==len(removed)==record['removed']
    a=np.r_[x[~delete],added];canonical=np.unique(a,axis=0)
    assert len(a)==len(canonical)==record['total']
    assert hashlib.sha256(canonical.astype('<i2').tobytes()).hexdigest()==record['canonical_point_sha256']
    np.savez_compressed(output,rational=a,radical=np.zeros_like(a),radicand=np.int32(2),squared_norm=np.int32(32),inner_product_limit=np.int32(16))
    print(json.dumps(dict(dimension=d,points=len(a),output=str(output),canonical_point_sha256=record['canonical_point_sha256'])))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('dimension',type=int);p.add_argument('source',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.dimension,a.source,a.output)
