#!/usr/bin/env python3
"""Build the enlarged exact 17D ground set, including both norm-six types.

Coordinates are (v[0:16], sqrt(2)*v[16])/6, with squared norm 8.
The 4320-point odd BW equator and two poles are fixed, not search vertices.
"""
import itertools
import json
import os
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA = HERE / 'data'


def span(basis):
    out = [0]
    for b in basis:
        out += [w ^ int(b) for w in out]
    return sorted(out)


def prepare():
    from ortools.sat.python import cp_model
    DATA.mkdir(exist_ok=True)
    octads = json.loads((ROOT / 'verifications/improved/dim18-bent-hexagon/data/octads.json').read_text())['octads']
    masks = [sum(1 << i for i in s) for s in octads]
    rm1 = sorted(set([0, 65535] + masks))
    assert len(rm1) == 32
    rm2 = [w for w in range(65536) if all((w & o).bit_count() % 2 == 0 for o in masks)]
    assert len(rm2) == 2048
    supports = [w for w in rm2 if w.bit_count() == 6]
    assert len(supports) == 448
    vectors, keys = [], []
    for support in supports:
        positions = [i for i in range(16) if support >> i & 1]
        for bits in range(64):
            word = sum(((bits >> j) & 1) << i for j, i in enumerate(positions))
            v = [0] * 17
            for i in positions:
                v[i] = -6 if word >> i & 1 else 6
            v[16] = 6
            vectors.append(v)
            keys.append((0, support, word))
    for q in range(16):
        for word in rm2:
            v = [-3 if word >> i & 1 else 3 for i in range(16)] + [6]
            v[q] *= 3
            vectors.append(v)
            keys.append((1, q, word))
    for word in rm2:
        vectors.append([-4 if word >> i & 1 else 4 for i in range(16)] + [4])
        keys.append((2, 0, word))
    V = np.asarray(vectors, dtype=np.int16)
    keys = np.asarray(keys, dtype=np.int32)
    assert V.shape == (63488, 17)
    assert np.all((V[:, :16].astype(np.int32)**2).sum(1) + 2*V[:, 16].astype(np.int32)**2 == 288)
    # Match the octads' labelling, rather than importing a separately labelled biplane.
    m0 = supports[0]
    family = set(w ^ m0 for w in rm1 if (w ^ m0).bit_count() == 6)
    assert len(family) == 16
    assert all((a & b).bit_count() == 2 for a, b in itertools.combinations(family, 2))
    A = [i for i, (kind, support, word) in enumerate(keys.tolist())
         if kind == 0 and support in family and word.bit_count() % 2 == 1]
    K = [w for w in rm2 if (w & m0).bit_count() % 2 == 0]
    unseen = set(K)
    fibers = []
    while unseen:
        r = min(unseen)
        fiber = sorted(r ^ h for h in rm1)
        fibers.append(fiber)
        unseen.difference_update(fiber)
    model = cp_model.CpModel()
    x = [model.new_bool_var(f'f{i}') for i in range(len(fibers))]
    for i, j in itertools.combinations(range(len(fibers)), 2):
        if any((fibers[i][0] ^ b).bit_count() == 4 for b in fibers[j]):
            model.add(x[i] + x[j] <= 1)
    model.maximize(sum(x))
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 1
    assert solver.solve(model) == cp_model.OPTIMAL
    B = sorted(w for i, f in enumerate(fibers) if solver.value(x[i]) for w in f)
    assert len(A) == 512 and len(B) == 192
    Bset = set(B)
    shift = next(s for s in K if Bset.isdisjoint(w ^ s for w in B))
    index = {tuple(k): i for i, k in enumerate(keys.tolist())}
    bp = [index[(2, 0, w)] for w in B]
    bm = [index[(2, 0, w ^ shift)] for w in B]
    full = np.concatenate([V, V * np.asarray([1]*16 + [-1], dtype=np.int16)])
    seed = np.asarray(A + bp + [len(V) + i for i in A + bm], dtype=np.int32)
    S = full[seed].astype(np.int32)
    gram = S[:, :16] @ S[:, :16].T + 2 * np.outer(S[:, 16], S[:, 16])
    np.fill_diagonal(gram, -1)
    assert gram.max() == 144
    # Independently check the entire candidate pool against the fixed equator.
    E = []
    for a, b in itertools.combinations(range(16), 2):
        for sa, sb in itertools.product((-12, 12), repeat=2):
            v = [0]*16
            v[a], v[b] = sa, sb
            E.append(v)
    for octad in octads:
        for bits in range(256):
            if bits.bit_count() % 2:
                v = [0]*16
                for j, p in enumerate(octad):
                    v[p] = -6 if bits >> j & 1 else 6
                E.append(v)
    E = np.asarray(E, dtype=np.int32)
    assert len(E) == 4320
    worst = max(int((V[a:a+512, :16].astype(np.int32) @ E.T).max()) for a in range(0, len(V), 512))
    assert worst == 144
    # Both poles: sqrt(2)*12/6, and |axis coefficient| <= 6.
    assert np.abs(full[:, 16]).max() == 6
    with (DATA / 'vertices17.bin').open('wb') as f:
        np.asarray(full.shape, dtype=np.uint32).tofile(f)
        full.tofile(f)
    np.savez_compressed(DATA / 'pool17.npz', vertices=full, keys=keys, seed=seed,
                        rm1=rm1, rm2=rm2, equator=E, family=sorted(family))
    (DATA / 'seed17.txt').write_text('\n'.join(map(str, seed)) + '\n')
    report = dict(vertices=len(full), half_vertices=len(V), fixed_points=4322,
                  seed_size=len(seed), seed_total=4322+len(seed),
                  candidate_equator_max_scaled_dot=worst,
                  seed_max_scaled_dot=int(gram.max()),
                  types_per_half=dict(six=28672, dual=32768, flat=2048))
    (DATA / 'preparation17.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    prepare()
