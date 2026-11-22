"""Verify rational heads, recount each affine shell, and enumerate short differences."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import json
from math import gcd, isqrt
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

HERE = Path(__file__).resolve().parent
CORE = HERE.parent / "dim40-affine-centre-refill"
sys.path.insert(0, str(CORE))
from verify_expanded import AFFINE, enclosing_form, load_module


def norm_head(head):
    return sum(int(x)**2 for x in head)


def head_dot(left, right):
    return sum(int(x)*int(y) for x, y in zip(left, right))


def body_cutoff(r, head, scale, s, other, denominator):
    """Upper bound on twice the integer body product, after moving the head term."""
    n = r*scale**2 + norm_head(head)
    m = s*denominator**2 + norm_head(other)
    return (isqrt(9*n*m) - 6*head_dot(head, other)) // (scale*denominator)


def head_rank(rows):
    pivots = {}
    for row in rows:
        v = [int(x) for x in row]
        for j in sorted(pivots):
            a = v[j]
            if a:
                u = pivots[j]
                v = [u[j]*x-a*y for x, y in zip(v, u)]
        if any(v):
            factor = gcd(*v)
            pivots[next(j for j, x in enumerate(v) if x)] = [x//factor for x in v]
    return len(pivots)


def verify(dimension, workers, compiler):
    assert __debug__, "Run with Python assertions enabled"
    with np.load(HERE / f"configuration_{dimension}.npz", allow_pickle=False) as data:
        assert all(data[key].dtype.kind in "iu" for key in data)
        g, residues, radii, heads, weights, points, tiny_heads, pure = (
            data[k].astype(np.int64) for k in ("gram", "residues", "radii", "heads",
            "populations", "tiny_bodies", "tiny_heads", "pure"))
        d, e = int(data["denominator"]), int(data["tiny_denominator"])
    assert g.shape == (40, 40) and np.array_equal(g, g.T) and np.abs(g).max() < 2**20
    assert not np.any(g.diagonal() % 2) and load_module(AFFINE/"verify.py").positive_definite(g) == 6561
    width = heads.shape[1]
    assert d > 0 and e > 0 and heads.shape == (len(radii), width)
    assert head_rank(np.vstack((heads, tiny_heads, pure))) == dimension-40
    assert residues.shape == (len(radii), 40) and np.all((0 <= residues) & (residues < 3))
    assert not np.any((residues @ g) % 3)
    assert weights.shape == radii.shape and set(radii) <= {10, 12, 14, 16, 18}
    assert points.shape == (len(tiny_heads), 40) and tiny_heads.shape[1] == pure.shape[1] == width
    assert np.abs(points).max(initial=0) < 128 and np.all((points @ g * points).sum(1) == 30)
    assert not np.any((points @ g) % 3)
    assert len(set(zip(map(tuple, points), map(tuple, tiny_heads)))) == len(points)
    assert all(norm_head(p) > 0 for p in pure)

    assert np.all(weights > 0)

    # Congruence modulo 18 follows from evenness and Gc = 0 mod 3.
    # Orient short differences up to sign and share their enumerations across all pairs.
    jobs = set()
    empty = set()
    refinements = defaultdict(list)
    groups = defaultdict(list)
    for i, (x, h) in enumerate(zip(points, tiny_heads)):
        groups[(tuple(x % 3), tuple(h))].append(i)

    def difference_keys(left, right, bound):
        delta = (left - right) % 3
        signed = tuple(map(int, delta))
        neg = tuple(map(int, -delta % 3))
        key = min(signed, neg)
        start = int(delta @ g @ delta) % 18 or 18
        return [((key, n), 1 if key == signed else -1) for n in range(start, bound, 18)]

    counts = {}
    for i, (c, radius, z, w) in enumerate(zip(residues, radii, heads, weights)):
        r = int(radius)
        key = (min(tuple(c), tuple(-c % 3)), 3*r)
        jobs.add(key)
        assert key not in counts or counts[key] == int(w)
        counts[key] = int(w)
        for j in range(i+1):
            s = int(radii[j])
            bound = 3*(r+s) - body_cutoff(r,z,d,s,heads[j],d)
            if i != j and r == s and np.array_equal(c, residues[j]):
                # Every body occurs at both heads, including pairs with zero difference.
                assert bound <= 0
            for key, _ in difference_keys(c, residues[j], bound):
                empty.add(key)
        for (cl, h), indices in groups.items():
            bound = 3*(r+10) - body_cutoff(r,z,d,10,h,e)
            if r == 10 and np.array_equal(c, cl):
                assert bound <= 0
            for key, sign in difference_keys(c, np.asarray(cl), bound):
                refinements[key].append((indices, sign, r))

    jobs.update(empty)
    jobs.update(refinements)
    jobs = sorted(jobs, key=lambda job: -job[1])
    print(f"Recounting {len(counts)} affine shells and testing {len(empty | refinements.keys())} short difference shells.", flush=True)
    with tempfile.TemporaryDirectory(prefix="kn_affine_heads_") as tmp:
        work = Path(tmp)
        binary = work/"enumerate"
        subprocess.run([compiler, "-O3", "-std=c++17", str(CORE/"enumerate_core.cpp"),
                        "-o", str(binary)], check=True)
        inputs = {}
        forms = {}
        enclosure = (HERE / "p48p_enclosure.npz" if dimension == 43 else
                     HERE.parent / "dim40-p48q-affine-refill/enclosure.npz")
        for norm in {job[1] for job in jobs}:
            scale, diag, lower, bound = enclosing_form(g, enclosure, norm)
            forms[norm] = scale, diag, lower, bound
            inputs[norm] = f"{scale} {bound}\n" + "\n".join(
                " ".join(map(str, row)) for row in (diag, *lower, *g)) + "\n"

        def prefixes(key):
            c, norm = key
            if counts.get(key, 0) < 100_000 or workers == 1:
                yield ()
                return
            scale, diag, lower, bound = forms[norm]

            def visit(chosen, remaining):
                if len(chosen) == 3:
                    yield chosen
                    return
                j = 39 - len(chosen)
                off = sum(int(lower[39-i,j])*value for i, value in enumerate(chosen))
                radius = isqrt(remaining // int(diag[j]))
                lo = -((radius+off)//scale)
                hi = (radius-off)//scale
                for value in range(lo+(int(c[j])-lo)%3, hi+1, 3):
                    cost = int(diag[j])*(scale*value+off)**2
                    yield from visit((*chosen, value), remaining-cost)
            yield from visit((), bound)

        def one(task):
            i, (key, prefix) = task
            c, n = key
            path, out = work/f"{i}.in", work/f"{i}.i8"
            path.write_text(inputs[n] + " ".join(map(str, c)) + "\n", encoding="ascii")
            cmd = [str(binary), str(path), str(out), "--norm", str(n), *map(str, prefix)]
            count = int(subprocess.check_output(cmd))
            assert out.stat().st_size == 40*count
            if key in empty:
                assert count == 0, (key, count)
            uses = refinements.get(key, ())
            if uses and count:
                v = np.fromfile(out, dtype=np.int8).reshape(-1, 40)
                for ids, sign, r in uses:
                    # v is oriented as family point minus explicit point.
                    dots = points[ids] @ g @ v.T
                    assert not np.any(2*sign*dots == 3*(r-10)-n)
            path.unlink()
            out.unlink()
            return key, count

        totals = Counter()
        tasks = ((key, prefix) for key in jobs for prefix in prefixes(key))
        with ThreadPoolExecutor(max_workers=workers) as pool:
            for key, count in pool.map(one, enumerate(tasks)):
                totals[key] += count
        assert all(totals[key] == expected for key, expected in counts.items())
        refined_vectors = sum(totals[key] for key in refinements)

    # Explicit points have the same body radius. Compute their body products
    # exactly in int64, and their rational-head thresholds in Python integers.
    grouped = defaultdict(list)
    for i, h in enumerate(tiny_heads):
        grouped[tuple(h)].append(i)
    items = list(grouped.items())
    weighted = points @ g
    for i, (h, ids) in enumerate(items):
        for j, (other, js) in enumerate(items[:i+1]):
            products = 2*(weighted[ids] @ points[js].T)
            if i == j:
                np.fill_diagonal(products, -100)
            assert np.all(products <= body_cutoff(10,h,e,10,other,e))
    for i, p in enumerate(pure):
        for q in pure[:i]:
            x = 2*head_dot(p,q)
            assert x <= 0 or x*x <= norm_head(p)*norm_head(q)
        for r,z,den in [*((int(r),z,d) for r,z in zip(radii,heads)), *((10,h,e) for h in grouped)]:
            x = 2*head_dot(p,z)
            assert x <= 0 or x*x <= norm_head(p)*(r*den**2+norm_head(z))

    total = int(weights.sum()) + len(points) + len(pure)
    return dict(status="PASS", dimension=dimension, points=total, complete_families=len(radii),
                recounted_body_points=int(weights.sum()), explicit_tiny=len(points), pure_heads=len(pure),
                distinct_shell_enumerations=len(jobs), difference_vectors_in_refinements=refined_vectors,
                explicit_tiny_pairs=len(points)*(len(points)-1)//2)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("dimensions", nargs="*", type=int, help="dimensions to replay (default: 41, 42, 43)")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--compiler", default="g++")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    expected = {41: 1_460_066, 42: 1_991_475, 43: 2_921_073}
    for dimension in args.dimensions or expected:
        if dimension not in expected:
            parser.error(f"unsupported dimension {dimension}")
        result = verify(dimension, args.workers, args.compiler)
        assert result["points"] == expected[dimension]
        print(json.dumps(result, indent=2), flush=True)
