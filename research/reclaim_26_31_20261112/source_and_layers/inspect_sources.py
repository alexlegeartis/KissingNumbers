"""Hash, schema, norm, layer, and 26D quadratic-block inspection."""
import argparse
import json
from pathlib import Path
import numpy as np
from source_paths import PACKINGSTAR_DEFAULT, SOURCE_AUDIT, packingstar_file
from source_union import greater3

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('--packingstar', type=Path, default=PACKINGSTAR_DEFAULT)
ap.add_argument('--output', type=Path)
args = ap.parse_args()
rows = []
for d in range(26, 32):
    path = packingstar_file(args.packingstar, d)
    with np.load(path, allow_pickle=False) as f:
        key = 'rational_points' if d == 26 else 'points'
        den = int(f['rational_denominator' if d == 26 else 'coordinate_denominator'])
        x = f[key].astype(np.int64)
        assert np.all(np.sum(x*x, axis=1) == 32*den*den)
        values, counts = np.unique(np.sum(x[:, :24]**2, axis=1)//(den*den), return_counts=True)
        row = dict(dimension=d, rational_points=len(x), rational_core_norm_counts=dict(zip(map(str, values), map(int, counts))))
        if d == 26:
            a, b = f['quadratic_A'].astype(np.int64), f['quadratic_B'].astype(np.int64)
            assert den == 1 and np.all(f['quadratic_denominators'] == 1) and np.all(f['quadratic_radicands'] == 3)
            assert np.all((a*a+3*b*b).sum(1) == 32) and np.all((a*b).sum(1) == 0)
            assert not greater3(x@a.T, x@b.T, 16).any()
            p, q = a@a.T+3*b@b.T, a@b.T+b@a.T
            np.fill_diagonal(p, -32)
            assert not greater3(p, q, 16).any()
            row.update(quadratic_points=len(a), exact_quadratic_cross_pairs=len(x)*len(a),
                       exact_quadratic_pairs=len(a)*(len(a)-1)//2)
        rows.append(row)
result = dict(packingstar_commit=SOURCE_AUDIT['PackingStar_commit'], dimensions=rows,
              scope='Exact hashes/norms/layers and 26D radical pairs; rational-rational angles are not checked here')
if args.output:
    args.output.write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
