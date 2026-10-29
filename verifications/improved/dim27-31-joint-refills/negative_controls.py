"""Require the geometric verifier to reject three distinct corruptions."""
from pathlib import Path
import tempfile
import numpy as np
from verify_geometry import verify


def main():
    a = np.zeros(27, dtype=np.int16)
    a[:8] = 4
    b = a.copy()
    b[0] *= -1
    assert int(a.astype(int) @ b.astype(int)) == 96 > 64
    cases = [('Distinct points', np.array([a, a])),
             ('Squared norms', np.array([a, a+1])),
             ('Unpaired-point angle', np.array([a, b]))]
    checks = 0
    with tempfile.TemporaryDirectory() as directory:
        for expected, points in cases:
            path = Path(directory)/'bad.npz'
            np.savez_compressed(path, points=points, squared_norm=128,
                                inner_product_limit=64, coordinate_denominator=2)
            for precision in (32, 64):
                try:
                    verify(path, precision)
                except ValueError as error:
                    assert str(error) == expected, str(error)
                    checks += 1
                else:
                    raise AssertionError('Corruption accepted: '+expected)
    print('PASS: %d rejection controls' % checks)


if __name__ == '__main__':
    main()
