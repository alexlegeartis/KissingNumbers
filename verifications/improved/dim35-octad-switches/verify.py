"""Replay the 415,358-point kissing configuration in dimension 35."""
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path

import numpy as np

P = Path(__file__).resolve().parent
SOURCE = P.parent / "dim35-shared-root-heads"
N = 35


def require_shape(X, magnitudes):
    assert X.ndim == 2 and X.shape[1] == N
    assert (np.sort(abs(X), axis=1) == np.array(magnitudes)).all()
    assert len(X) == len(set(map(tuple, X)))


def signs(template, centers):
    c, *leaves = template
    assert len(set(template)) == 5 and min(template) >= 0 and max(template) < N
    for a in centers:
        for v in product((-2, 2), repeat=4):
            if v.count(-2) % 2:
                continue
            x = [0] * N
            x[c] = a
            for j, q in zip(leaves, v):
                x[j] = q
            yield tuple(x)


def mask(indices):
    assert len(indices) == len(set(indices)) == 8
    assert all(type(j) is int and 0 <= j < N for j in indices)
    return sum(1 << j for j in indices)


def maxdot(A, B, cap, within=False):
    M = A @ B.T
    if within:
        np.fill_diagonal(M, 0)
    d = int(M.max())
    assert d <= cap, (d, cap)
    return d


def replay():
    if not __debug__:
        raise RuntimeError("Run without -O")
    c = json.loads((SOURCE / "code_data.json").read_text())
    w = json.loads((P / "witness.json").read_text())
    # Historical sparse points and their original root ledger are frozen inputs.
    f = SOURCE / "shared_positive.npz"
    assert sha256(f.read_bytes()).hexdigest() == "9d4b35ac5d18f40011282ec2e83753bdffee947351584a77a910baef875768cf"
    with np.load(f, allow_pickle=False) as source:
        old, heads, positions, old_removed = (source[key].astype(np.int64) for key in
            ("new_stars", "heads", "dense_positions", "removed_roots"))
    removal = w["remove"]
    a, b = removal["start"], removal["stop"]
    assert (a, b) == (62, 70)
    assert set(map(tuple, old[a:b])) == set(signs(removal["template"], (4,)))
    inserts = [x for t in w["both"] for x in signs(t, (-4, 4))]
    inserts += list(signs(w["positive"], (4,)))
    assert len(inserts) == 104
    stars = np.concatenate((old[:a], old[b:], np.array(inserts)), axis=0)
    require_shape(stars, [0]*30+[2]*4+[4])
    require_shape(heads, [1]*33+[3]*2)
    assert len(stars) == 246 and len(heads) == 12
    assert len(set(positions)) == len(positions) == 32 and 0 <= min(positions) <= max(positions) < N
    assert (abs(heads[:, positions]) == 1).all()

    original = [int(m) for m, _ in c["full_parity_supports"]]
    octads = original.copy()
    swapped = set()
    for i, before, after in w["octads"]:
        assert 0 <= i < len(octads) and i not in swapped and octads[i] == mask(before)
        octads[i] = mask(after)
        swapped.add(i)
    assert len(swapped) == 7 and len(set(octads)) == len(octads) == 2202
    assert len(set(octads) - set(original)) == 7
    assert all(0 < m < 1 << N and m.bit_count() == 8 for m in octads)
    assert all((A & B).bit_count() <= 4 for A, B in combinations(octads, 2))
    prior_parities = dict(json.loads((SOURCE / "full_parities.json").read_text()))
    assert set(prior_parities) == set(original) and set(prior_parities.values()) <= {0, 1}
    for i, m in enumerate(octads):
        ix = [j for j in range(N) if m >> j & 1]
        parity = 0 if i in swapped else prior_parities[m]
        assert (2 * abs(stars[:, ix]).sum(axis=1) <= 16).all()
        h = heads[:, ix]
        best = 2*abs(h).sum(axis=1) - 4*abs(h).min(axis=1)*(((h < 0).sum(axis=1) % 2) != parity)
        assert (best <= 20).all()
    # Octads on one support have fixed parity, hence at least two sign flips.
    # Root/root and root/octad bounds are at most 16 by support and signs.
    allroots = []
    for i, j in combinations(range(N), 2):
        for a, b in product((-4, 4), repeat=2):
            x = [0] * N
            x[i], x[j] = a, b
            allroots.append(x)
    A = np.array(allroots, dtype=np.int64)
    blocked = (stars @ A.T > 16).any(axis=0) | (heads @ A.T > 20).any(axis=0)
    roots = A[~blocked]
    D = set(map(tuple, A[blocked]))
    old_D = set(map(tuple, old_removed))
    assert len(D) == 208 and len(old_D) == 134
    assert (len(D - old_D), len(old_D - D)) == (80, 6)
    assert len(roots) == 2172
    pair_maxima = dict(star_star=maxdot(stars, stars, 16, True),
                       head_head=maxdot(heads, heads, 25, True),
                       star_head=maxdot(stars, heads, 20),
                       root_star=maxdot(roots, stars, 16),
                       root_head=maxdot(roots, heads, 20))

    G = list(map(int, c["dense_generator"]))
    shift = int(c["dense_coset_shift"])
    assert len(G) == 17 and all(0 < g < 1 << 32 for g in G) and 0 <= shift < 1 << 32
    B = [0]
    for g in G:
        B += [x ^ g for x in B]
    assert len(B) == len(set(B)) == 131072 and min(x.bit_count() for x in B[1:]) >= 8
    # Dense norm-square 32; min binary distance 8 gives dot at most 16.
    dense = 1 - 2*((np.array(B, dtype=np.int64)[:, None] ^ shift) >> np.arange(32) & 1)
    pair_maxima["head_dense"] = maxdot(heads[:, positions], dense, 20)
    assert (abs(stars[:, positions]).sum(axis=1) <= 16).all()
    # Dense against octad or root is at most 16 or 8, respectively.
    # The head norm-square is 51, all other layers have norm-square 32,
    # and layers have distinct support sizes 35, 32, 8, 5 and 2.
    total = len(dense) + 128*len(octads) + len(roots) + len(stars) + len(heads)
    assert total == 415358
    return dict(status="PASS", dimension=N, points=total, baseline=415336,
                gain=22, dense=len(dense), octads=128*len(octads), roots=len(roots),
                stars=len(stars), heads=len(heads), removed_roots=len(D),
                roots_restored=6, roots_newly_blocked=80, pair_maxima=pair_maxima)


if __name__ == "__main__":
    print(json.dumps(replay(), sort_keys=True))
