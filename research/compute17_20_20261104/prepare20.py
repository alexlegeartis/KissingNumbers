#!/usr/bin/env python3
"""20D orbit model with *all* sign words and both octad sign parities.

The 760 roots are fixed. The full ground set contains 2^20 dense vectors and
130*256 signed octads. Quotienting by a subgroup of the shortened Golay code
keeps exact orbit weights and all conflicts, while allowing support sign classes
and dense vectors outside the usual 4096-word dual code to change together.
"""
import argparse
import itertools
import json
from pathlib import Path
import sys
import numpy as np
from orbits17 import DATA, write_graph

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'common'))
from golay import golay24


def basis_of(words):
    pivots = {}
    for value in words:
        value = int(value)
        for p, b in sorted(pivots.items(), reverse=True):
            if value >> p & 1:
                value ^= b
        if value:
            pivots[value.bit_length()-1] = value
    return list(pivots.values())


def reduction(basis, bits):
    pivots = {}
    for value in basis:
        value = int(value)
        for p, b in sorted(pivots.items(), reverse=True):
            if value >> p & 1:
                value ^= b
        if value:
            pivots[value.bit_length()-1] = value
    def reduce(values):
        values = np.asarray(values, dtype=np.uint32).copy()
        for p, b in sorted(pivots.items(), reverse=True):
            values ^= ((values >> p) & 1) * np.uint32(b)
        return values
    nonpivots = [p for p in range(bits) if p not in pivots]
    def coordinates(values):
        values = reduce(values)
        out = np.zeros_like(values)
        for i, p in enumerate(nonpivots):
            out |= ((values >> p) & 1) << i
        return out
    return reduce, coordinates, nonpivots


def make(rank, seed):
    G = golay24().astype(np.uint32)
    powers = np.uint32(1) << np.arange(20, dtype=np.uint32)
    punctured = G[:, 4:] @ powers
    shortened = punctured[G[:, :4].sum(1) == 0]
    assert len(shortened) == 256
    blocks = shortened[np.asarray([int(w).bit_count() == 8 for w in shortened])]
    assert len(blocks) == 130
    rng = np.random.default_rng(seed)
    basis = basis_of(rng.permutation(shortened))[:rank]
    assert len(basis) == rank
    H = [0]
    for b in basis:
        H += [w ^ b for w in H]
    assert min(w.bit_count() for w in H if w) >= 8
    reduce, coordinates, nonpivots = reduction(basis,20)
    ndense = 2**(20-rank)
    ids = np.arange(ndense, dtype=np.uint32)
    representatives = np.zeros(ndense, dtype=np.uint32)
    for i,p in enumerate(nonpivots):
        representatives |= ((ids >> i) & 1) << p
    assert np.array_equal(coordinates(representatives), ids)
    # All dense-dense conflicts: Hamming distance <=4, modulo H.
    errors = [sum(1 << p for p in positions) for k in range(1,5)
              for positions in itertools.combinations(range(20),k)]
    connections = np.unique(coordinates(errors))
    assert 0 not in connections
    rows = [list(map(int, i ^ connections)) for i in range(ndense)]
    weights = [2**rank]*ndense
    # The zero set in the first removed Golay coordinate is a 2048-word
    # distance-at-least-five code, and contains H, so it descends to this quotient.
    flat_seed = punctured[G[:,0] == 0]
    assert len(flat_seed) == 2048
    selected = list(map(int, np.unique(coordinates(flat_seed))))
    block_offsets, block_reps, block_weights = [], [], []
    for block in blocks:
        positions = [p for p in range(20) if int(block) >> p & 1]
        restricted = [sum(((h >> p) & 1) << i for i,p in enumerate(positions)) for h in H]
        rbasis = basis_of(restricted)
        rr, rc, rnonpivots = reduction(rbasis,8)
        local_patterns = np.arange(256, dtype=np.uint32)
        rids = rc(local_patterns)
        count = 2**(8-len(rbasis))
        assert set(map(int,rids)) == set(range(count))
        local_reps = np.asarray([int(local_patterns[rids == j].min()) for j in range(count)])
        offset = len(rows)
        block_offsets.append(offset)
        block_reps.append(local_reps)
        block_weights.append(2**len(rbasis))
        for pattern in local_reps:
            assert int(pattern).bit_count() % 2 == int(rr([pattern])[0]).bit_count() % 2
            # Flipping a single sign conflicts within the signed octad cube.
            adjacent = np.unique(rc(int(pattern) ^ (np.uint32(1) << np.arange(8, dtype=np.uint32))))
            assert rids[int(pattern)] not in adjacent
            rows.append([offset+int(k) for k in adjacent])
            weights.append(2**len(rbasis))
            if int(pattern).bit_count() % 2:
                selected.append(offset+int(rids[int(pattern)]))
        restriction = np.zeros(ndense, dtype=np.uint32)
        for j,p in enumerate(positions):
            restriction |= ((representatives >> p) & 1) << j
        # An octad vector conflicts with a dense vector precisely when all eight signs agree.
        neighbor_ids = offset + rc(restriction)
        for i,j in enumerate(neighbor_ids):
            rows[i].append(int(j))
            rows[int(j)].append(i)
    rows = [np.asarray(sorted(set(row)), dtype=np.uint32) for row in rows]
    name = f'q20_r{rank}_s{seed}'
    edges = write_graph(DATA / (name+'.graph'), rows)
    chosen = np.zeros(len(rows), dtype=bool)
    chosen[selected] = True
    assert all(not chosen[rows[i]].any() for i in selected)
    total = 760+sum(weights[i] for i in selected)
    assert total == 19448
    np.savez_compressed(DATA / (name+'.npz'), weights=weights, seed=selected, subgroup=H,
                        dense_representatives=representatives, blocks=blocks,
                        block_offsets=block_offsets, block_reps=np.asarray(block_reps,dtype=object),
                        block_weights=block_weights, fixed=760, dimension=20)
    report = dict(name=name, dimension=20, vertices=len(rows), directed_edges=edges,
                  dense_orbits=ndense, dense_connections=len(connections),
                  full_dense_candidates=1048576, full_octad_candidates=33280,
                  group_rank=rank, orbit_weights=sorted(set(weights)), baseline=total,
                  subgroup_seed=seed)
    (DATA / (name+'.json')).write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--rank',type=int,default=8)
    p.add_argument('--seed',type=int,default=1)
    a=p.parse_args()
    make(a.rank,a.seed)
