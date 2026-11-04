"""Pin public inputs and put all seeds in convenient, documented coordinates."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data'


def normalize(x):
    return x / np.linalg.norm(x, axis=1, keepdims=True)


def stats(x):
    x = normalize(x)
    g = x @ x.T
    vals = g[np.triu_indices(len(x), 1)]
    v = np.maximum(vals - .5, 0)
    return dict(n=len(x), d=x.shape[1], max_cos=float(vals.max()),
                sum_squared_violation=float(v @ v),
                violating_pairs=int(np.count_nonzero(v)),
                pairs_above_05001=int(np.count_nonzero(vals > .5001)),
                min_squared_distance=float(2 - 2*vals.max()))


def main():
    manifest = {}
    for p in sorted(DATA.glob('*.txt')):
        manifest[p.name] = dict(bytes=p.stat().st_size,
                               sha256=hashlib.sha256(p.read_bytes()).hexdigest())
    canonical = normalize(np.loadtxt(DATA / 'canonical840.txt'))
    record = normalize(np.loadtxt(DATA / 'gram841_coordinates.txt'))
    u = np.loadtxt(DATA / 'alignment_U.txt')
    perm = np.loadtxt(DATA / 'alignment_permutation.txt', dtype=int) - 1
    assert np.array_equal(np.sort(perm), np.arange(840))
    # Paper's original_841_index and pi_of_i columns are 1-based.
    old = np.delete(record, 750, axis=0) @ u
    aligned = np.empty((841, 12))
    aligned[perm] = old
    aligned[-1] = record[750] @ u
    err = np.sum((aligned[:840] - canonical)**2, axis=1)
    assert abs(err.sum() - 11.96241654955465) < 1e-6
    assert err[:816].sum() < 4.667
    np.save(DATA / 'canonical840.npy', canonical)
    np.save(DATA / 'record841.npy', normalize(aligned))
    baseline = dict(canonical840=stats(canonical), record841=stats(aligned),
                    alignment_squared_error=float(err.sum()))
    for n in [842, 844]:
        gram = np.loadtxt(DATA / f'Gram{n}.txt')
        assert gram.shape == (n, n)
        assert np.max(abs(gram-gram.T)) < 1e-10
        w, v = eigh(gram, subset_by_index=[n-12, n-1])
        x = normalize(v * np.sqrt(w))
        np.save(DATA / f'published{n}.npy', x)
        baseline[f'published{n}'] = stats(x)
        baseline[f'published{n}']['gram_reconstruction_error'] = float(np.max(abs(x@x.T-gram)))
    (ROOT / 'inputs.json').write_text(json.dumps(dict(
        sources=dict(quaternionic='https://github.com/k-nic/quaternionic',
                     near_misses='https://github.com/k-nic/841_in_12D/tree/main/beyond841',
                     paper='https://arxiv.org/abs/2609.09179'),
        files=manifest, baseline=baseline), indent=2)+'\n')
    print(json.dumps(baseline, indent=2))


if __name__ == '__main__':
    main()
