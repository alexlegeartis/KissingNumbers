"""Small adversarial inputs for the complete coordinate checker."""
from pathlib import Path
import json
import tempfile
import numpy as np
from verify_geometry import verify


def main():
    a = np.zeros(30, dtype=np.int16)
    a[:8] = 4
    b = a.copy()
    b[0] = -4  # Same squared norm 128, but forbidden product 96.
    good = np.asarray([a, -a])
    cases = {
        'wrong_norm': (np.asarray([a, np.zeros(30, dtype=np.int16)]), 128),
        'duplicate': (np.asarray([a, a]), 128),
        'same_norm_bad_angle': (np.asarray([a, b]), 128),
        'wrong_metadata': (good, 32),
    }
    rejected = []
    with tempfile.TemporaryDirectory() as directory:
        for name, (points, norm) in cases.items():
            path = Path(directory) / (name + '.npz')
            np.savez_compressed(path, points=points, squared_norm=norm,
                                inner_product_limit=64, coordinate_denominator=2)
            try:
                verify(path)
            except ValueError:
                rejected.append(name)
            else:
                raise ValueError('Malformed input was accepted: ' + name)
        path = Path(directory) / 'valid_control.npz'
        np.savez_compressed(path, points=good, squared_norm=128,
                            inner_product_limit=64, coordinate_denominator=2)
        if verify(path)['status'] != 'PASS':
            raise ValueError('Valid control rejected')
    return dict(status='PASS', rejected=rejected, valid_control='PASS')


if __name__ == '__main__':
    print(json.dumps(main(), indent=2))
