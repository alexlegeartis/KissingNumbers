"""Exact Delsarte–Goethals dense-owner replay in dimension 65."""
from itertools import combinations
from functools import reduce
from pathlib import Path
from fractions import Fraction
from collections import Counter
import hashlib
import importlib.util
import json

HERE = Path(__file__).resolve().parent
V = HERE.parent
PARENT = V / "dim64-65-signed-parity" / "construction.py"
COUNTS = V / "dim64-65-partial-nonary-signs" / "claims.json"

def parent_module():
    spec = importlib.util.spec_from_file_location("owner_parent", PARENT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def gf32mul(a, b):
    q = 0
    while b:
        if b & 1:
            q ^= a
        a <<= 1
        if a & 32:
            a ^= 0x25  # X^5 + X^2 + 1
        b >>= 1
    return q

def tr(a):
    t = 0
    for _ in range(5):
        t ^= a
        a = gf32mul(a, a)
    assert t in (0, 1)
    return t

BASIS = list(combinations(range(6), 2))

def form(a, b):
    # Symmetric binary bilinear form on F32 and the alternating extension to six bits.
    def bil(x, y):
        return tr(gf32mul(a, gf32mul(x, y))
                  ^ gf32mul(b, gf32mul(gf32mul(x, x), y) ^ gf32mul(x, gf32mul(y, y))))
    B = [[bil(1 << i, 1 << j) for j in range(5)] for i in range(5)]
    return sum((B[i][j] ^ (B[i][i] & B[j][j])) << k if j < 5 else B[i][i] << k
               for k, (i, j) in enumerate(BASIS))

def span(rows):
    X = [0]
    for r in rows:
        X += [r ^ w for w in X]
    return X

def main():
    assert __debug__, "Assertions must be enabled"
    assert hashlib.sha256(PARENT.read_bytes()).hexdigest() == "3af72f18fb50eefb71d044b7624aad31897d2edcd50e4955993873a83d1ca9ec"
    assert hashlib.sha256(COUNTS.read_bytes()).hexdigest() == "73d6a1f6b7e5b06bc0885748b80bd02781ab85d833b599b2afa82fd4e4d92db0"
    P = parent_module()
    parent, _ = P.extended_bch(13)
    old_owners, _ = P.extended_bch(23)
    parent_checks = [sum(((a & b).bit_count() % 2) << i for i, b in enumerate(parent))
                     for a in parent]
    dense = [reduce(int.__xor__, (r for j, r in enumerate(parent) if w >> j & 1), 0)
             for w in P.kernel(parent_checks, len(parent))]
    D, piv = P.rref(dense, 64)
    assert len(D) == 28 and len(P.rref(old_owners, 64)[1]) == 16
    assert all(P.in_span(x, D, piv) for x in old_owners)

    X = [1]
    for _ in range(1, 63):
        X.append(P.gf64_multiply(X[-1], 2))
    X.append(0)
    assert sorted(X) == list(range(64))
    one = (1 << 64) - 1
    A = [sum(((x >> i) & 1) << j for j, x in enumerate(X)) for i in range(6)]
    Q = [A[i] & A[j] for i, j in BASIS]
    assert len(P.rref([one, *A, *Q], 64)[1]) == 22
    assert all(P.in_span(w, D, piv) for w in [one, *A, *Q])
    assert all(gf32mul(x, 1) == x and (x == 0 or len({gf32mul(x, y) for y in range(32)}) == 32)
               for x in range(32))
    forms = sorted({form(a, b) for a in range(32) for b in range(32)})
    assert len(forms) == 1024
    diff = {a ^ b for a, b in combinations(forms, 2)}
    affine = span([one, *A])
    assert len(affine) == 128 and one in affine
    evalquad = lambda m: reduce(int.__xor__, (v for j, v in enumerate(Q) if m >> j & 1), 0)
    mindist = min((evalquad(m) ^ w).bit_count() for m in diff for w in affine)
    assert mindist == 24 and min(w.bit_count() for w in affine if w) == 32
    new_count = len(forms) * len(affine)
    assert new_count == 131072

    # Sparse owner conditions in the inherited proof: 336 affine GF(4)-line supports,
    # eight even-parity signs each. Distinct owners in one cap have raw dot <= 16.
    lines = {tuple(sorted([x, x ^ d, x ^ P.gf4_times_alpha(d),
                                  x ^ d ^ P.gf4_times_alpha(d)]))
             for x in range(64) for d in range(1, 64)}
    pairs = [pair for line in lines for pair in combinations(line, 2)]
    assert len(lines) == 336 and len(pairs) == len(set(pairs)) == 2016
    sparse = len(lines) * 8
    assert 64 - 2*mindist == 4*4 == 16
    assert Fraction(3, 4)*16 + 16 < 32
    assert Fraction(3, 4)*64 - 16 == 32 and 8*4 == 32

    values = {r["dimension"]: r["lower_bound"] for r in json.loads(COUNTS.read_text())["dimensions"]}
    assert values == {64: 335795354, 65: 335863580}
    assert values[65] == values[64] + (1 << len(old_owners)) + sparse + 2
    total = values[64] + new_count + sparse + 2
    result = {"dimension": 65, "previous": values[65], "lower_bound": total,
              "gain": total-values[65], "new_dense_owners": new_count,
              "old_dense_owners": 1 << len(old_owners), "sparse_owners": sparse,
              "quadratic_difference_masks": len(diff), "minimum_hamming_distance": mindist}
    assert result == json.loads((HERE / "claims.json").read_text())
    return result

if __name__ == "__main__":
    print(json.dumps(main(), indent=2))
