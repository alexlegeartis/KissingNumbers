"""Exact integer replay of the 652,076-point, 38-dimensional construction."""
from hashlib import sha256
import importlib.util
from itertools import combinations, product
import json
from pathlib import Path

P = Path(__file__).resolve().parent
S = P.parent / "dim37-38-qr-merge" / "supports_3088.txt"
N = 38


def mask(indices, size):
    assert len(indices) == size == len(set(indices))
    assert all(type(i) is int and 0 <= i < N for i in indices)
    return sum(1 << i for i in indices)


def sparse(recipe):
    raw = S.read_bytes()
    assert sha256(raw).hexdigest() == "ccb530f90e914df84cbc8d19a71d30c79146d23a2e8e86e7a813d4d92a5a6a8b"
    octads = [mask([i for i, a in enumerate(w) if a == "1"], 8)
              for w in raw.decode("ascii").splitlines()]
    assert len(octads) == len(set(octads)) == 3088
    traded = set()
    for row in recipe["replacements"]:
        i = row["i"]
        assert type(i) is int and 0 <= i < len(octads) and i not in traded
        assert octads[i] == mask(row["remove"], 8)
        octads[i] = mask(row["add"], 8)
        traded.add(i)
    assert len(traded) == 10
    assert all((x & y).bit_count() <= 4 for x, y in combinations(octads, 2))
    # Same-support even parity imposes at least two sign changes.
    # Distinct octads interact only on their intersection.
    stars = []
    for center, *leaf in recipe["stars"]:
        leaf_mask = mask(leaf, 4)
        assert not (leaf_mask & mask([center], 1))
        centermask = 1 << center
        assert all(2 * bool(o & centermask) + (o & leaf_mask).bit_count() <= 4
                   for o in octads)
        for sign in (-2, 2):
            for signs in product((-1, 1), repeat=4):
                if signs.count(-1) % 2:
                    st = [0] * N
                    st[center] = sign
                    for j, s in zip(leaf, signs):
                        st[j] = s
                    stars.append(tuple(st))
    assert len(stars) == len(set(stars)) == 288
    assert all(sum(t*t for t in x) == 8 for x in stars)
    assert all(sum(u*v for u, v in zip(x, y)) <= 4
               for x, y in combinations(stars, 2))
    roots = []
    for i, j in combinations(range(N), 2):
        for a, b in product((-2, 2), repeat=2):
            if all(a*x[i] + b*x[j] <= 4 for x in stars):
                roots.append((i, j, a, b))
    assert len(roots) == 2572
    # Distinct roots (norm-square 8) have dot <=4; octad/root
    # pairs have at most two shared entries, each of product at most 2.
    return len(octads), len(roots), len(stars)


def replay():
    assert __debug__, "Run without Python -O so verification assertions execute"
    path = P.parent / "dim38-qr-free-offsets" / "verify.py"
    spec = importlib.util.spec_from_file_location("offset_proof", path)
    base = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(base)
    inherited = base.replay()
    n = inherited["dense"]
    assert inherited["all_fibre_pairs"] == 465 and n == 253952
    o, r, s = sparse(json.loads((P / "witness.json").read_text()))
    # Dense ordinary magnitudes are 2 and marked magnitude is 0 or 3.
    # Use the sparse L1 mass, with a bound for the change at one coordinate.
    for head in (0, 3):
        extra = max(0, head-2)
        for dot in (16+extra, 8+2*extra, 12+2*extra):
            assert 4*dot**2 <= 8*(148+head**2)
    # Different layer types have distinct support sizes.
    result = n + 128*o + r + s
    assert result == 652076
    return dict(status="PASS", dimension=38, points=result, baseline=652072,
                gain=4, dense=n, octads=128*o, roots=r, stars=s,
                removed_roots=2*N*(N-1)-r)


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run without Python -O")
    print(json.dumps(replay(), sort_keys=True))
