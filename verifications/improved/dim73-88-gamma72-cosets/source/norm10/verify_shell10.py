"""Independent integer verification of the new 73D norm-ten shell layer."""
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def main():
    start = perf_counter()
    parent=HERE.parent/'dim73_surgery'
    if not parent.exists(): parent=ROOT/'research/nonlattice73_20261025'
    lattice_path=parent/'lattice/certificate.npz'
    if not lattice_path.exists(): lattice_path=ROOT/'research/lattice73_20261025/result/certificate.npz'
    lattice = np.load(lattice_path)
    g = lattice['gamma_gram'].astype(np.int64)
    x = lattice['heads'].astype(np.int64)
    u = np.load(parent/'surgery_class.npy').astype(np.int64)
    path = HERE/'shell10_selected.npy'
    y = np.load(path).astype(np.int64)
    assert y.shape[1] == 72 and len(y)>0 and len(np.unique(y,axis=0)) == len(y)
    yg = y@g
    assert np.all(np.sum(yg*y,axis=1) == 10)
    largest = -10
    smallest = 10
    headmax = capmax = 0
    for i in range(0,len(y),128):
        dots = yg[i:i+128]@y.T
        dots[np.arange(len(dots)),i+np.arange(len(dots))] = 0
        largest = max(largest,int(dots.max()))
        smallest = min(smallest,int(dots.min()))
        headmax = max(headmax,int(np.abs(yg[i:i+128]@x.T).max()))
        capmax = max(capmax,int(np.abs(yg[i:i+128]@u.T).max()))
    assert largest<=3 and smallest>=-6 and headmax<=5 and capmax<=3
    # Norm: (16/25)*10+8/5 = 8. For distinct oriented heads the
    # upper-point dot lies between -56/25 and 88/25, so antipodal completion
    # is also safe. Equator products <=(4/5)*5=4 by minimum 8 of Gamma72.
    assert -100 <= 16*smallest+40 <= 16*largest+40 <= 100
    assert 2*headmax+10 <= 20  # (2/5)*D + 2 <= 4.
    # Mixed old caps: 6*sqrt(3)+4*sqrt(5)<20, checked by two squarings.
    assert 400-108-80 > 0 and 48**2*15 < (400-108-80)**2
    previous = 6218175600+2*len(x)+2*len(u)
    report = {'status':'PASS','dimension':73,'oriented_heads':len(y),'added_points':2*len(y),
              'previous_bound':previous,'bound':previous+2*len(y),
              'minimum_distinct_head_product':smallest,'maximum_distinct_head_product':largest,
              'maximum_coset_product':headmax,'maximum_cap_product':capmax,
              'integer_pair_products':len(y)*(len(y)+len(x)+len(u)),
              'certificate_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
              'seconds':perf_counter()-start,
              'inherited_inputs':'Previously verified 73D non-lattice certificate and published Gamma72 minimum/shell size.'}
    (HERE/'shell10_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':
    main()
