#!/usr/bin/env python3
"""Exact subgroup quotients of the enlarged 17D graph, without assuming full RM(1,4)."""
import argparse
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE / 'data'


def read_graph(path):
    head = np.fromfile(path, dtype=np.uint32, count=2)
    n, magic = map(int, head)
    assert magic == 0x4b4e4731
    offsets = np.memmap(path, dtype=np.uint64, mode='r', offset=8, shape=(n+1,))
    neighbors = np.memmap(path, dtype=np.uint32, mode='r', offset=8+8*(n+1), shape=(int(offsets[-1]),))
    return offsets, neighbors


def write_graph(path, rows):
    offsets = np.zeros(len(rows)+1, dtype=np.uint64)
    offsets[1:] = np.cumsum([len(row) for row in rows], dtype=np.uint64)
    with Path(path).open('wb') as f:
        np.asarray([len(rows), 0x4b4e4731], dtype=np.uint32).tofile(f)
        offsets.tofile(f)
        for row in rows:
            np.asarray(row, dtype=np.uint32).tofile(f)
    return int(offsets[-1])


def make(rank, constant):
    z = np.load(DATA / 'pool17.npz')
    rm1 = list(map(int, z['rm1']))
    basis = []
    generated = {0}
    for w in ([65535] if constant else []) + [w for w in rm1 if w not in (0, 65535)]:
        if w not in generated:
            # For the nonconstant case, explicitly avoid generating the constant.
            enlarged = generated | {w ^ a for a in generated}
            if not constant and 65535 in enlarged:
                continue
            basis.append(w)
            generated = enlarged
            if len(basis) == rank:
                break
    assert len(basis) == rank and ((65535 in generated) == constant)
    H = sorted(generated)
    keys = z['keys'].tolist()
    index = {tuple(k): i for i, k in enumerate(keys)}
    half = len(keys)
    orbit = np.full(2*half, -1, dtype=np.int32)
    representatives, sizes = [], []
    for i in range(2*half):
        if orbit[i] >= 0:
            continue
        layer, local = divmod(i, half)
        kind, support, word = keys[local]
        members = sorted({index[(kind, support, word ^ (h & support if kind == 0 else h))] + layer*half for h in H})
        oid = len(representatives)
        assert all(orbit[j] < 0 for j in members)
        orbit[members] = oid
        representatives.append(i)
        sizes.append(len(members))
    assert min(sizes) == max(sizes) == 2**rank
    off, nb = read_graph(DATA / 'graph17.bin')
    rows = []
    for oid, rep in enumerate(representatives):
        row = np.unique(orbit[nb[off[rep]:off[rep+1]]])
        assert oid not in row, 'an orbit contains an internal conflict'
        rows.append(row)
    name = f'q17_r{rank}_c{int(constant)}'
    edges = write_graph(DATA / (name+'.graph'), rows)
    seed = np.unique(orbit[z['seed']])
    assert sum(sizes[i] for i in seed) == len(z['seed'])
    np.savez_compressed(DATA / (name+'.npz'), orbit=orbit, representatives=representatives,
                        weights=sizes, seed=seed, subgroup=H, fixed=4322, dimension=17)
    report = dict(name=name, vertices=len(rows), directed_edges=edges,
                  orbit_size=2**rank, baseline=5730, subgroup=H)
    (DATA / (name+'.json')).write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--rank', type=int, default=4)
    p.add_argument('--constant', type=int, choices=(0,1), default=1)
    a = p.parse_args()
    make(a.rank, bool(a.constant))
