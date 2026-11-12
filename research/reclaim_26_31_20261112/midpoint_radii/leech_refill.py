"""Complete free Leech-core insertion scan for an explicit changed seed.

This checks a finite candidate family, not the full geometry of the seed.
The Leech shell is regenerated from the repository's Golay construction.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'source_and_layers'))
from source_paths import fresh_leech_shell


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('source', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    z = np.load(args.source, allow_pickle=False)
    x = z['points'].astype(np.int16)
    denominator = int(z['coordinate_denominator'])
    if not np.all(np.sum(x.astype(np.int64)**2, axis=1) == 32*denominator**2):
        raise ValueError('Seed norm and coordinate scale')
    leech = fresh_leech_shell()
    on_core = np.all(x[:, 24:] == 0, axis=1)
    if np.any(x[on_core, :24] % denominator):
        raise ValueError('Nonintegral equator')
    present = {row.tobytes() for row in (x[on_core, :24]//denominator)}
    full = {row.tobytes() for row in leech}
    if not present <= full:
        raise ValueError('Equator not contained in the generated Leech shell')
    missing = leech[[row.tobytes() not in present for row in leech]]
    # Zero-tail candidates see only the distinct old core projections.
    heads = np.unique(x[~on_core, :24], axis=0)
    if 24*int(np.abs(leech).max())*int(np.abs(heads).max()) >= 2**24:
        raise ValueError('Exact binary32 accumulation guard')
    alive = np.ones(len(missing), dtype=bool)
    products = 0
    for start in range(0, len(heads), 512):
        ids = np.flatnonzero(alive)
        if not len(ids):
            break
        block = heads[start:start+512]
        dots = missing[ids].astype(np.float32) @ block.astype(np.float32).T
        products += len(ids)*len(block)
        alive[ids] = np.all(dots <= 16*denominator, axis=1)
    additions = np.zeros((int(alive.sum()), x.shape[1]), dtype=np.int16)
    additions[:, :24] = missing[alive]*denominator
    # Independently replay every survivor against every source point.
    for point in additions:
        dots = np.sum(x.astype(np.int64)*point.astype(np.int64), axis=1)
        if np.any(dots > 16*denominator**2):
            raise ValueError('Direct integer survivor check')
    np.savez_compressed(args.out/'additions.npz', points=additions)
    report = {
        'source': str(args.source.resolve()),
        'source_sha256': hashlib.sha256(args.source.read_bytes()).hexdigest(),
        'baseline': len(x), 'missing_leech_points': len(missing),
        'distinct_nonzero_tail_core_projections': len(heads),
        'added': len(additions), 'total': len(x)+len(additions),
        'products': products, 'seconds': time.monotonic()-started,
        'scope': 'Complete free insertion from the fixed Leech core only; no seed-validity or exchange optimality claim.'
    }
    (args.out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
