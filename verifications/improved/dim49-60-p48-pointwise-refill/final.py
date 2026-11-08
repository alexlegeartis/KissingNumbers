"""Derive registered counts from the supplied final ownership certificates."""
from pathlib import Path
import json
import numpy as np

D = Path(__file__).resolve().parent
POLES = (2, 6, 12, 24, 40, 72, 126, 240, 306, 510, 604, 841)
CLASSES = (1, 3, 6, 12, 20, 36, 63, 120, 153, 255, 302, 420)


def results():
    metadata = json.loads((D / 'certificate.json').read_text())
    rows = metadata['rows']
    if [r['dimension'] for r in rows] != list(range(49, 61)):
        raise ValueError('Expected one ordered certificate per dimension 49 through 60')
    out = []
    for row, poles, classes in zip(rows, POLES, CLASSES):
        with np.load(D / ('heads%d.npz' % row['dimension']), allow_pickle=False) as data:
            bits = data['retained_bits']
            extra = data['extra_points']
            offsets = data['extra_offsets']
            if bits.dtype != np.uint8 or bits.shape != (classes, 973):
                raise ValueError('Invalid retained-mask type or shape')
            masks = np.unpackbits(bits, axis=1)
            if masks[:, 7783:].any() or extra.ndim != 2 or extra.shape[1] != 48:
                raise ValueError('Invalid padding or extra-point shape')
            if offsets.shape != (classes + 1,) or offsets[0] != 0 or offsets[-1] != len(extra) or np.any(np.diff(offsets) < 0):
                raise ValueError('Invalid extra-point partition')
            lines = int(masks[:, :7783].sum()) + len(extra)
        total = 52416000 + 2 * lines + poles
        if (lines != row['selected_lines'] or total != row['lower_bound']
                or poles != row['poles'] or classes != row['classes']
                or total <= row['previous_repository_bound']):
            raise ValueError('Derived count disagrees with certificate metadata')
        out.append(dict(dimension=row['dimension'], total=total))
    return out


if __name__ == '__main__':
    for row in results():
        print(json.dumps(row))
