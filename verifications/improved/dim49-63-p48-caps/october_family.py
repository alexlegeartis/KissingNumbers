"""Compare complete independently disjoint P48 families, without mixing them."""
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent


def load_refilled():
    """Expand exact generator words and retained-index bitmaps into coordinates."""
    with np.load(ROOT / 'october/family_refilled.npz', allow_pickle=False) as archive:
        data = {key: archive[key] for key in archive.files}
    assert all(value.dtype.kind in 'iu' for value in data.values())
    references, words = data['transform_references'], data['transform_words']
    bits, base = data['retained_bits'], data['base']
    assert references.shape == (577,) and words.shape == (577, 30)
    assert bits.dtype == np.uint8 and bits.shape == (577, (len(base)+7)//8)
    masks = np.unpackbits(bits, axis=1)
    assert not masks[:, len(base):].any()
    masks = masks[:, :len(base)]
    indices = [np.flatnonzero(row).astype(np.int32) for row in masks]
    sizes = np.array([len(row) for row in indices], dtype=np.int64)
    assert np.array_equal(sizes, data['sizes'])
    with np.load(ROOT / 'october/family_55_61.npz', allow_pickle=False) as archive:
        historical = archive['transforms'].copy()
    generators = np.load(ROOT / 'data/p48p_gens.npy').astype(np.int64)
    transforms = []
    for reference, word in zip(references, words):
        if reference >= 0:
            assert reference < len(historical) and np.all(word == -1)
            transform = historical[reference].copy()
        else:
            assert reference == -1 and np.all((word >= 0) & (word < len(generators)))
            transform = np.eye(48, dtype=np.int64)
            for index in word:
                generator = generators[index]
                assert 48*int(np.abs(transform).max())*int(np.abs(generator).max()) < 2**63
                transform = transform @ generator
        transforms.append(transform)
    return dict(base=base, transforms=np.array(transforms),
                retained=np.concatenate(indices), offsets=np.r_[0, np.cumsum(sizes)],
                projection=data['projection'], sizes=sizes, gram=data['gram'])


def families():
    old = np.load(ROOT / 'data/class_sizes_p48.npy').ravel()
    with np.load(ROOT / 'october/family_55_61.npz', allow_pickle=False) as data:
        new = data['sizes'].copy()
    with np.load(ROOT / 'october/family_refilled.npz', allow_pickle=False) as data:
        refilled = data['sizes'].copy()
    return {'class_sizes_p48.npy': np.sort(old)[::-1],
            'october/family_55_61.npz': np.sort(new)[::-1],
            'october/family_refilled.npz': np.sort(refilled)[::-1]}


def best_prefix(count):
    options = [(int(sizes[:count].sum()), name) for name, sizes in families().items()
               if len(sizes) >= count]
    assert options
    return max(options)
