"""Verify the parent lattice, the three empty projected shells, and the 44D configuration."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from functools import lru_cache
import json
from math import isqrt, prod
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
import sympy as sp

from moments import HERE, K, verify_families

CORE = HERE.parent / "dim40-affine-centre-refill"
sys.path.insert(0, str(CORE))
from verify_expanded import AFFINE, enclosing_form, load_module


def verify_minimum(gram, workers, compiler):
    scale, diagonal, lower, bound = enclosing_form(gram, HERE/"ambient.npz", 4)
    dimension = len(gram)
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        binary, data = work/"enumerate", work/"input"
        subprocess.run([compiler, "-O3", "-std=c++17", "-DKN_SHELL_RANK=48",
                        "-DKN_RESIDUE_MODULUS=1", str(CORE/"enumerate_core.cpp"),
                        "-o", str(binary)], check=True)
        with data.open("w") as f:
            print(scale, bound, file=f)
            for row in (diagonal, *lower, *gram, [0]*dimension):
                print(*row, file=f)

        def prefixes(prefix, remaining):
            if len(prefix) == 4:
                yield prefix
                return
            k = dimension-1-len(prefix)
            offset = sum(int(lower[dimension-1-i,k])*x for i,x in enumerate(prefix))
            limit = isqrt(remaining//int(diagonal[k]))
            lo = -((limit+offset)//scale)
            if not any(prefix):
                lo = max(lo, 0)  # A nonzero vector and its negative have the same norm.
            for x in range(lo, (limit-offset)//scale+1):
                cost = int(diagonal[k])*(scale*x+offset)**2
                yield from prefixes((*prefix,x), remaining-cost)

        def count(prefix):
            cmd = [str(binary), str(data), "/dev/null", "--through", "4", *map(str,prefix)]
            return int(subprocess.check_output(cmd))

        with ThreadPoolExecutor(max_workers=workers) as pool:
            assert sum(pool.map(count, prefixes((),bound))) == 0


def pair_section_size():
    # For (a,b)=3 on shell six, P selects (x,a)=3. On that slice x and
    # a-x both have norm six, so (x,b) belongs to {-3,0,1,2,3,6}.
    x = sp.Symbol("x")
    p = sp.Poly(prod(Fraction(1,3-i)*(x-i) for i in range(-3,3)), x)
    q = sp.Poly(prod(Fraction(1,-i)*(x-i) for i in (-3,1,2,3,6)), x)

    @lru_cache(None)
    def pairings(n,m):
        if n == 0 and m == 0:
            return 1
        if n == 0:
            return 6*(m-1)*pairings(0,m-2)
        within = 6*(n-1)*pairings(n-2,m) if n >= 2 else 0
        across = 3*m*pairings(n-1,m-1) if m else 0
        return within+across

    result = 0
    for (n,), c in p.terms():
        for (m,), e in q.terms():
            if (n+m) % 2 == 0:
                result += c*e*Fraction(52_416_000 * 6**((n+m)//2) * pairings(n,m),
                                       prod(range(48,48+n+m,2)))
    # At x=±a the scalar product with b is ±3, both roots of Q.
    assert result == 1656
    return int(result)


def verify(workers, compiler):
    assert __debug__, "Run with Python assertions enabled"
    with np.load(HERE/"ambient.npz", allow_pickle=False) as data:
        assert all(data[k].dtype.kind in "iu" for k in data)
        gram, normals, points = (data[k].astype(np.int64) for k in ("gram","normals","section_points"))
    assert gram.shape == (48,48) and np.array_equal(gram,gram.T)
    assert not np.any(gram.diagonal() % 2) and abs(gram).max() < 2**16
    assert load_module(AFFINE/"verify.py").positive_definite(gram) == 1
    assert normals.shape == (4,48) and abs(normals).max() < 128
    assert np.array_equal(normals @ gram @ normals.T, K)
    verify_minimum(gram, workers, compiler)
    calibration, forms = verify_families()
    assert points.shape == (pair_section_size(),48) and abs(points).max() < 128
    assert len(np.unique(points,axis=0)) == len(points)
    products = points @ gram
    assert np.all((products*points).sum(axis=1) == 6)
    profiles = products @ normals.T
    assert np.all(profiles[:,0] == 3) and np.all(profiles[:,1] == 0)
    # Each profile has first entries (3,0) or (-3,0), so these counts
    # follow from the certified slice and its image under negation.
    counts = []
    for row in calibration:
        assert abs(row[0]) == 3 and row[1] == 0
        oriented = row if row[0] == 3 else -row
        counts.append(int(np.all(profiles == oriented, axis=1).sum()))
    sizes = forms @ [1,*counts]
    with np.load(HERE/"families.npz", allow_pickle=False) as data:
        assert np.array_equal(sizes, data["populations"]) and np.all(sizes > 0)
    total = int(sizes.sum())
    assert counts == [0,0,0,134] and total == 4_940_961
    return dict(status="PASS", dimension=44, points=total, complete_families=len(sizes),
                parent_minimum=6, section_points=len(points), calibration_counts=counts)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--compiler", default="g++")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    print(json.dumps(verify(args.workers, args.compiler), indent=2))
