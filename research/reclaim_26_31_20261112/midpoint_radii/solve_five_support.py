"""Replay the shared-deletion MILP used for the six-thousand-point pilot.

The search result is numerical. Every proposed positive exchange receives an
independent direct-integer changed-pair check before a proposal file is written.
This does not replace full geometric certification of any eventual record.
"""
from pathlib import Path
import argparse
import importlib.util
import json
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=150)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[3]
    path = root/'research/exceptional_tail_replacement_20261108/optimize_rational.py'
    spec = importlib.util.spec_from_file_location('_five_support_optimizer', path)
    optimizer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(optimizer)
    z = np.load(args.directory/'pool.npz', allow_pickle=False)
    candidates = z['candidates'].astype(np.int16)
    ptr, old = z['blocker_indptr'], z['blocker_indices']
    if not np.all(np.sum(candidates.astype(np.int64)**2,axis=1)==32):
        raise ValueError('Candidate norm')
    blockers = [old[ptr[i]:ptr[i+1]] for i in range(len(candidates))]
    report = optimizer.optimize(candidates, blockers, norm=32, seconds=args.seconds)
    if report['gain'] > 0:
        source = np.load(args.source, allow_pickle=False)
        points = source['points'].astype(np.int16)
        if int(source['squared_norm']) == 128:
            if np.any(points % 2):
                raise ValueError('Source not integral at norm 32')
            points //= 2
        if not np.all(np.sum(points.astype(np.int64)**2,axis=1)==32):
            raise ValueError('Source norm')
        added = candidates[report['added']]
        keep = np.ones(len(points), dtype=bool)
        keep[report['removed']] = False
        final = np.concatenate([points[keep], added])
        if len(np.unique(final, axis=0)) != len(final):
            raise ValueError('Duplicate final points')
        for index, point in enumerate(added):
            dots = np.sum(final.astype(np.int64)*point.astype(np.int64),axis=1)
            dots[int(keep.sum())+index] = -32
            if np.any(dots > 16):
                raise ValueError('Changed-pair inequality')
        np.savez_compressed(args.directory/'proposal.npz', points=final,
                            squared_norm=32, inner_product_limit=16,
                            coordinate_denominator=1)
        report['independent_changed_pairs'] = 'PASS'
    (args.directory/'optimization.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
