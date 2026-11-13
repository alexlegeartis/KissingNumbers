"""Bounded exact E8-octad reflections of the November 31D source.

The complete image is isometric to the source. Selection is restricted by
the explicit blocker cutoff; negative outcomes are not global bounds.
"""
from pathlib import Path
from itertools import product
import argparse
import hashlib
import json
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'verifications/improved/dim27-31-joint-refills'))
from search_sextets import refill, save


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('source', type=Path)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--trials', type=int, default=8)
    ap.add_argument('--cutoff', type=int, default=32)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    z = np.load(args.source, allow_pickle=False)
    x = z['points'].astype(np.int16)
    norm = int(z['squared_norm'])
    if norm == 32:
        x *= 2
    elif norm != 128:
        raise ValueError('Source normalization')
    normal = [0, 2, 4, 5, 6, 10, 11, 23]
    roots = [np.asarray(s, dtype=np.int16) for s in product((-1, 1), repeat=8)
             if s[0] == 1 and np.prod(s) == 1]
    rng = np.random.default_rng(20261113)
    order = [63] + [i for i in rng.permutation(len(roots)) if i != 63]
    baseline = x.copy()
    started = time.monotonic()
    records = []
    for trial, rid in enumerate(order[:args.trials]):
        root = roots[rid]
        scalar = x[:, normal].astype(np.int32) @ root.astype(np.int32)
        if np.any(scalar % 4):
            raise ValueError('This exact reflection requires a larger denominator')
        y = x.copy()
        y[:, normal] -= ((scalar // 4)[:, None] * root).astype(np.int16)
        if not np.all(np.sum(y.astype(np.int64)**2, axis=1) == 128):
            raise ValueError('Reflected norm')
        x, rec = refill(x, y, args.cutoff)
        rec.update(trial=trial, root=root.tolist(), total=len(x), seconds=time.monotonic()-started)
        records.append(rec)
        print(json.dumps(rec), flush=True)
        if len(x) > len(baseline):
            save(args.out / f'coordinates31_{len(x)}.npz', x)
            old = {tuple(row) for row in baseline}
            new = {tuple(row) for row in x}
            np.savez_compressed(args.out / 'exchange.npz',
                                added=np.asarray(sorted(new-old), dtype=np.int16),
                                removed=np.asarray(sorted(old-new), dtype=np.int16))
        report = dict(source=str(args.source), source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
                      baseline=len(baseline), final=len(x), cutoff=args.cutoff, records=records,
                      scope='Finite whole-isometry unions; full output replay required for any claimed bound.')
        (args.out / 'report.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
