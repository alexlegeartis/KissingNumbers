"""Count the supplied nearest-subset, repaired and scaled Gamma72 lifts."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def results():
    rows = {}
    folder = ROOT / 'source/dim61_96'
    for name in ['gamma72_coset_certificate.json','gamma72_repaired_lifts.json',
                 'gamma72_scaled_4.json','gamma72_scaled_8.json','gamma72_scaled_16.json']:
        path = folder / name
        if not path.exists():
            continue
        report = json.loads(path.read_text(encoding='utf-8'))
        assert name == 'gamma72_coset_certificate.json' or report['status'].startswith('Verified:')
        for dim, row in report['dimensions'].items():
            if 'retained_head_lines' in row:
                count = row['retained_head_lines']
                count = sum(count) if isinstance(count, list) else count
                total = 6218175600 + 4*count
            else:
                total = 6218175600 + 2*row['head_subset_size']*row['cosets'] + row['pure_tail_points']
            assert total == row['bound']
            rows[int(dim)] = max(rows.get(int(dim), 0), total)
    import numpy as np
    surgery = ROOT / 'source/dim73_surgery'
    with np.load(surgery / 'lattice/certificate.npz') as certificate:
        nearest_count = len(certificate['heads'])
    class_count = len(np.load(surgery / 'surgery_class.npy'))
    rows[73] = max(rows.get(73, 0), 6218175600 + 2*nearest_count + 2*class_count)
    return [dict(dimension=dim, total=total) for dim,total in sorted(rows.items())]


if __name__ == '__main__':
    for row in results():
        print(json.dumps(row))
