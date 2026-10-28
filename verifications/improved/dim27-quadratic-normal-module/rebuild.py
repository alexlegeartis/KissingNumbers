"""Rebuild the full certificate from pinned PackingStar data and our module."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent


def main(source,output):
    meta=json.loads((ROOT/'data/provenance.json').read_text())
    assert hashlib.sha256(source.read_bytes()).hexdigest()==meta['packingstar_file_sha256']
    with np.load(source,allow_pickle=False) as z:
        den=int(z['coordinate_denominator']);assert np.all(z['points']%den==0)
        x=(z['points']//den).astype(np.int16)
    assert x.shape==(209586,27)
    complement=np.setdiff1d(np.arange(24),meta['normal_octad'])
    sector=np.all(x[:,complement]==0,axis=1);assert int(sector.sum())==370
    with np.load(ROOT/'data/normal_module.npz',allow_pickle=False) as z:
        ma=z['rational'];mb=z['radical'];assert ma.shape==mb.shape==(384,27)
    a=np.r_[x[~sector],ma];b=np.r_[np.zeros_like(x[~sector]),mb]
    canonical=np.unique(np.c_[a,b].astype('<i2'),axis=0)
    assert len(canonical)==209600
    assert hashlib.sha256(canonical.tobytes()).hexdigest()==meta['canonical_coefficient_sha256']
    np.savez_compressed(output,rational=a,radical=b,radicand=np.int32(2),squared_norm=np.int32(32),inner_product_limit=np.int32(16))
    print(json.dumps(dict(points=len(a),output=str(output),canonical_coefficient_sha256=meta['canonical_coefficient_sha256'])))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('--output',type=Path,default=Path('rebuilt27.npz'));a=p.parse_args();main(a.source,a.output)
