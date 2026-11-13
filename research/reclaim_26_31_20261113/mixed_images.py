"""Exact mixed core/tail images of the pinned November 31D configuration."""
from pathlib import Path
import argparse
import hashlib
import importlib.util
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
    ap.add_argument('--trials', type=int, default=12)
    ap.add_argument('--cutoff', type=int, default=16)
    ap.add_argument('--seed', type=int, default=20261104)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    spec = importlib.util.spec_from_file_location('old_image_maps', ROOT / 'research/algebraic_switches_20261104/search_images.py')
    maps = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(maps)
    z = np.load(args.source, allow_pickle=False)
    x = z['points'].astype(np.int16)
    norm = int(z['squared_norm'])
    if norm == 32:
        x *= 2
    elif norm != 128:
        raise ValueError('Source normalization')
    baseline = x.copy()
    started = time.monotonic()
    records = []
    for trial, action in enumerate(maps.maps(x.shape[1], args.trials, args.seed, extended=True)):
        m = np.asarray(action['numerator'], dtype=np.int64)
        denominator = action['denominator']
        numerator = x.astype(np.int64) @ m
        if np.any(numerator % denominator):
            raise ValueError('Exact fractional image requires larger scaling')
        y = (numerator // denominator).astype(np.int16)
        if not np.all(np.sum(y.astype(np.int64)**2, axis=1) == 128):
            raise ValueError('Reflected norm')
        x, rec = refill(x, y, args.cutoff)
        rec.update(trial=trial, action=action, total=len(x), seconds=time.monotonic()-started)
        records.append(rec)
        print(json.dumps({k:v for k,v in rec.items() if k != 'action'}), flush=True)
        if len(x) > len(baseline):
            save(args.out / f'coordinates31_{len(x)}.npz', x)
            old, new = {tuple(row) for row in baseline}, {tuple(row) for row in x}
            np.savez_compressed(args.out / 'exchange.npz',
                                added=np.asarray(sorted(new-old), dtype=np.int16),
                                removed=np.asarray(sorted(old-new), dtype=np.int16))
        report = dict(source=str(args.source), source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
                      baseline=len(baseline), final=len(x), cutoff=args.cutoff, records=records,
                      scope='Finite whole-isometry unions; full output replay required for any claimed bound.')
        (args.out / 'report.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
