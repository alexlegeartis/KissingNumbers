"""Verify the shipped compact class family and bind it to the counted sizes."""
from pathlib import Path
import argparse
import json
import os
import runpy
import sys
import time
import numpy as np

HERE = Path(__file__).resolve().parent
if not __debug__:
    raise RuntimeError('Run the verifier with assertions enabled (without -O).')
sys.path.insert(0, str(HERE / 'gpu'))
sys.path.insert(0, str(HERE))
from cpu_certificate import verify

DATA = HERE.parent / 'data'
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--complete', action='store_true', help='also check all cap directions, partitions, poles and counts for 87-94D')
args = parser.parse_args()
started = time.perf_counter()
certificate = DATA / 'symmetry_classes.npz'
with np.load(certificate, allow_pickle=False) as z:
    sizes = z['sizes']
    assert np.array_equal(sizes, np.load(DATA / 'class_sizes_symmetry.npy'))
verify(certificate)
if args.complete:
    assert len(sizes) >= 16632, 'certificate does not cover all eight dimensions'
    for name in ('verify_parts.py', 'verify_poles.py', 'verify_record.py'):
        original_argv, original_cwd = sys.argv, Path.cwd()
        try:
            sys.argv = [str(HERE / name), '--k', '15', '16', '17', '18', '19', '20', '21', '22']
            os.chdir(HERE.parent)
            try:
                runpy.run_path(str(HERE / name), run_name='__main__')
            except SystemExit as result:
                if result.code not in (None, 0):
                    raise
        finally:
            sys.argv = original_argv
            os.chdir(original_cwd)
    rows = json.loads((HERE.parent / 'rows81-95.json').read_text())
    covered = [row for row in rows if 87 <= row['dim'] <= 94]
    assert sorted(row['dim'] for row in covered) == list(range(87, 95))
    sums = np.r_[0, np.cumsum(np.sort(sizes)[::-1].astype(np.int64))]
    for row in covered:
        assert row['k'] == row['dim'] - 72 and row['W'] == 'lattice'
        with np.load(DATA / ('parts_%d.npz' % row['k'])) as partition:
            assert row['T'] == len(partition['triples'])
            assert row['P'] == len(partition['pairs']) == 0
        assert row['T'] <= len(sizes)
        expected = 6218175600 + 4 * int(sums[row['T']]) + row['poles']
        assert row['value'] == expected, 'driver value is not this verified construction'
        print('VERIFIED tau(%d) >= %d' % (row['dim'], expected), flush=True)
    print('ALL EIGHT CAP CONSTRUCTIONS VERIFIED in %.3f CPU seconds' % (time.perf_counter()-started))
