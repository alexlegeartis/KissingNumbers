"""Reject wrong norms, duplicates and illegal signed or opposite-pair angles."""
from pathlib import Path
import tempfile
import numpy as np
from verify_geometry import verify


def save(path, points):
    np.savez_compressed(path, points=points, squared_norm=128,
                        inner_product_limit=64, coordinate_denominator=2)


if __name__ == '__main__':
    a = np.zeros(28, dtype=np.int16); a[:2] = 8
    b = np.zeros(28, dtype=np.int16); b[0] = 8; b[1:5] = 4
    bad_norm = a.copy(); bad_norm[0] += 1
    cases = {'norm': [a, bad_norm], 'duplicate': [a, a],
             'signed_angle': [a, b], 'opposite_pair_angle': [a, -a, b]}
    with tempfile.TemporaryDirectory() as directory:
        for name, rows in cases.items():
            path = Path(directory)/(name+'.npz'); save(path, np.asarray(rows))
            try:
                verify(path)
            except ValueError:
                print('PASS: rejected '+name)
            else:
                raise RuntimeError('Failed to reject '+name)
