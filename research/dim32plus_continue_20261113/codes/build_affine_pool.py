"""Recover AG(2,8) from the dual octads and generate its support orbit.

This is discovery preprocessing. Every retained result must be checked against
the actual binary generator and support/sign inequalities before registration.
"""
from pathlib import Path
import itertools
import json
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = ROOT / 'verifications/improved/dim61-68-95-96-code-layers/source'


def recover_plane():
    lines = list(map(int, json.loads((SOURCE/'shared_parity/dual_octads.json').read_text())))
    assert len(lines) == len(set(lines)) == 72
    remaining = set(lines)
    classes = []
    while remaining:
        a = min(remaining)
        group = sorted(h for h in remaining if h == a or not h & a)
        assert len(group) == 8 and sum(group) == (1 << 64)-1
        assert all(not a & b for a,b in itertools.combinations(group,2))
        classes.append(group)
        remaining.difference_update(group)
    assert len(classes) == 9
    rows, cols = classes[:2]
    grid = [[(a&b).bit_length()-1 for b in cols] for a in rows]
    assert all((a&b).bit_count() == 1 for a in rows for b in cols)
    origin = grid[0][0]
    diagonals = classes[2]
    through = next(h for h in diagonals if h >> origin & 1)
    p0 = [next(j for j,c in enumerate(cols) if through & a & c) for a in rows]
    inv = [p0.index(j) for j in range(8)]
    perms = [tuple(inv[next(j for j,c in enumerate(cols) if h & a & c)] for a in rows)
             for h in diagonals]
    identity = tuple(range(8))
    generated = {identity:0}
    for p in perms:
        if p in generated:
            continue
        bit = len(generated)
        assert bit in (1,2,4)
        add = {tuple(p[q[i]] for i in range(8)): value|bit for q,value in generated.items()}
        assert not set(add) & set(generated)
        generated.update(add)
    assert set(generated) == set(perms) and len(generated) == 8
    labels = [next(value for p,value in generated.items() if p[0] == i) for i in range(8)]
    point = np.zeros((8,8),dtype=np.int64)
    xy = np.zeros((64,2),dtype=np.int64)
    for i in range(8):
        for j in range(8):
            x,y = labels[i],labels[inv[j]]
            point[x,y] = grid[i][j]
            xy[grid[i][j]] = (x,y)
    mult = np.zeros((8,8),dtype=np.int64)
    for family in classes[1:]:
        h = next(h for h in family if h >> origin & 1)
        f = [next(y for y in range(8) if h >> int(point[x,y]) & 1) for x in range(8)]
        assert f[0] == 0
        assert all(f[a^b] == f[a]^f[b] for a in range(8) for b in range(8))
        mult[f[1]] = f
    assert mult[1].tolist() == list(range(8))
    assert np.array_equal(mult,mult.T)
    assert all(mult[a,mult[b,c]] == mult[mult[a,b],c] for a,b,c in itertools.product(range(8),repeat=3))
    assert all(mult[a,b^c] == mult[a,b]^mult[a,c] for a,b,c in itertools.product(range(8),repeat=3))
    assert all(mult[a,b] for a,b in itertools.product(range(1,8),repeat=2))
    return lines,classes,point,xy,mult


def main():
    lines,classes,point,xy,mult = recover_plane()
    with np.load(SOURCE/'shared_parity/global.npz',allow_pickle=False) as z:
        seed = int(z['new_supports'][0])
    positions = [i for i in range(64) if seed >> i & 1]
    sx,sy = xy[positions].T
    tx,ty = np.indices((8,8)).reshape(2,64)
    unit = np.left_shift(np.uint64(1),point.astype(np.uint64))
    orbit = set()
    matrices = 0
    for a,b,c,d in itertools.product(range(8),repeat=4):
        if mult[a,d] == mult[b,c]:
            continue
        matrices += 1
        x = mult[a,sx] ^ mult[b,sy]
        y = mult[c,sx] ^ mult[d,sy]
        masks = np.bitwise_or.reduce(unit[x[None,:]^tx[:,None],y[None,:]^ty[:,None]],axis=1)
        orbit.update(map(int,masks))
    assert matrices == 3528 and seed in orbit
    values = np.array(sorted(orbit),dtype=np.uint64)
    assert all(int(x).bit_count() == 16 for x in values)
    with np.load(SOURCE/'dim61_96/code64_supports.npz',allow_pickle=False) as z:
        old = z['words'].copy()
        generator = z['generator'].copy()
    checks = [sum(int(v)<<i for i,v in enumerate(row)) for row in generator]
    assert all(not((int(x)&h).bit_count()%2) for x in values for h in checks)
    values.tofile(HERE/'candidate_supports.bin')
    old.tofile(HERE/'old_supports.bin')
    report = dict(status='PASS',dual_octads=len(lines),parallel_classes=len(classes),
                  field_multiplication=mult.tolist(),point_coordinates=xy.tolist(),
                  seed=seed,linear_maps=matrices,affine_images=matrices*64,
                  distinct_candidates=len(values),old_supports=len(old),
                  all_candidates_in_parent_dual=True)
    (HERE/'affine_plane.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('field_multiplication','point_coordinates')}),flush=True)


if __name__ == '__main__':
    main()
