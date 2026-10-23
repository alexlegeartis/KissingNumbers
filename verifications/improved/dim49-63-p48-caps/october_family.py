"""Compare the two independently disjoint P48 class families, without mixing them."""
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent


def families():
    old = np.load(ROOT / 'data/class_sizes_p48.npy').ravel()
    with np.load(ROOT / 'october/family_55_61.npz', allow_pickle=False) as data:
        new = data['sizes'].copy()
    return {'class_sizes_p48.npy': np.sort(old)[::-1],
            'october/family_55_61.npz': np.sort(new)[::-1]}


def best_prefix(count):
    options = [(int(sizes[:count].sum()), name) for name, sizes in families().items()
               if len(sizes) >= count]
    assert options
    return max(options)
