"""Freshly enumerate the affine core and verify the explicit cross-layer pairs."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
import json
from math import isqrt
import os
from pathlib import Path
import subprocess
import tempfile

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")

import numpy as np

from verify import AFFINE, CORE_COUNT, GLUING, HERE, load_module, witness


def enclosing_form(gram, path, squared_norm):
    with np.load(path, allow_pickle=False) as enclosure:
        assert all(enclosure[key].dtype.kind in "iu" for key in enclosure.files)
        scale = int(enclosure["P"])
        diagonal = enclosure["D"].astype(object)
        lower = enclosure["M"].astype(object)
        inverse = enclosure["inverse_times_three"].astype(object)
    rank = len(gram)
    assert gram.shape == (rank, rank)
    assert 0 < scale < 2**40 and diagonal.shape == (rank,) and lower.shape == (rank, rank)
    assert inverse.shape == (rank, rank)
    assert np.array_equal(gram.astype(object) @ inverse, 3 * np.eye(rank, dtype=object))
    assert np.all(np.diag(lower) == scale) and not np.any(np.triu(lower, 1))
    assert np.all((0 < diagonal) & (diagonal < 2**36)) and np.max(np.abs(lower)) < 2**40
    error = scale**3 * gram.astype(object) - (lower * diagonal) @ lower.T
    eta = max(sum(abs(int(x)) for x in row) for row in error)
    trace = Fraction(sum(int(x) for x in np.diag(inverse)), 3)
    # If xQx <= squared_norm, then ||x||^2 <= squared_norm trace(Q^-1).
    # The residual of the scaled rational LDL form has operator norm <= eta.
    bound = squared_norm * (scale**3 + eta * trace)
    radius = -(-bound.numerator // bound.denominator)
    assert 0 < radius < 2**120
    return scale, diagonal, lower, radius


def check_minimum(gram, compiler, work, enclosure_path):
    scale, diagonal, lower, radius = enclosing_form(gram, enclosure_path, 4)
    data = work / "minimum.input"
    with data.open("w", encoding="ascii") as stream:
        print(40, scale, radius, file=stream)
        for row in (diagonal, *lower, *gram):
            print(*row, file=stream)
    binary = work / "exact_minimum"
    subprocess.run([compiler, "-O3", "-std=c++17", str(AFFINE / "exact_minimum.cpp"),
                    "-o", str(binary)], check=True)
    result = json.loads(subprocess.check_output([str(binary), str(data)], text=True))
    assert result["nonzero_vectors_of_norm_at_most_4"] == 0
    return result


def enumerate_core(gram, residue, compiler, work, *, enclosure_path=None, workers=1):
    path = enclosure_path if enclosure_path is not None else AFFINE / "enclosure.npz"
    scale, diagonal, lower, radius = enclosing_form(gram, path, 54)
    data = work / "core.input"
    with data.open("w", encoding="ascii") as stream:
        print(scale, radius, file=stream)
        for row in (diagonal, *lower, *gram, residue):
            print(*row, file=stream)
    binary = work / "enumerate_core"
    subprocess.run([compiler, "-O3", "-std=c++17", str(HERE / "enumerate_core.cpp"),
                    "-o", str(binary)], check=True)

    def prefixes(prefix, remaining):
        if len(prefix) == 3:
            yield prefix
            return
        k = 39 - len(prefix)
        shift = sum(int(lower[39 - i, k]) * value for i, value in enumerate(prefix))
        limit = isqrt(remaining // int(diagonal[k]))
        lo = -((limit + shift) // scale)
        hi = (limit - shift) // scale
        lo += (int(residue[k]) - lo) % 3
        for value in range(lo, hi + 1, 3):
            offset = scale * value + shift
            yield from prefixes((*prefix, value), remaining - int(diagonal[k]) * offset**2)

    def part(task):
        i, prefix = task
        result = work / f"core-{i}.i8"
        out = subprocess.check_output([str(binary), str(data), str(result), *map(str, prefix)],
                                      text=True)
        rows = np.fromfile(result, dtype=np.int8).reshape(-1, 40)
        assert int(out) == len(rows)
        return rows

    if workers == 1:
        core = part((0, ()))
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            core = np.concatenate(list(pool.map(part, enumerate(prefixes((), radius)))))
    assert len(core) == CORE_COUNT
    return core


def main():
    assert __debug__, "Run without Python -O"
    parser = argparse.ArgumentParser()
    parser.add_argument("--compiler", default="g++")
    args = parser.parse_args()
    certificate, auxiliary, refill, extra = witness()
    gram = np.array(json.loads((AFFINE / "gram.json").read_text()), dtype=np.int64)
    geometry = np.load(GLUING / "geometry.npz", allow_pickle=False)
    assert np.array_equal(geometry["body_gram"], gram)
    residue = geometry["body_residues"][certificate["centre"]] % 3
    parent = load_module(AFFINE / "verify.py")
    assert np.array_equal(gram, gram.T) and not np.any(np.diag(gram) % 2)
    assert parent.positive_definite(gram) == 6561
    with tempfile.TemporaryDirectory(prefix="kn40_affine_") as tmp:
        work = Path(tmp)
        parent.check_minimum(gram, args.compiler, work)
        core = enumerate_core(gram, residue, args.compiler, work)
    assert np.all(core % 3 == residue)
    points = [core, np.r_[auxiliary[0], extra], auxiliary[1], refill]
    norms = [54, 36, 30, 6]
    for a, norm in zip(points, norms):
        assert len(np.unique(a, axis=0)) == len(a)
        assert np.all((a @ gram * a).sum(axis=1) == norm)
    maxima = {}
    for i in range(len(points)):
        # The freshly checked rank-40 minimum proves the core/core bound.
        for j in range(max(1, i), len(points)):
            raw, _ = parent.max_product(points[i], points[j], gram, i == j, False)
            assert raw <= 0 or 4 * raw**2 <= norms[i] * norms[j]
            maxima[f"{norms[i]}/{norms[j]}"] = raw
    total = sum(map(len, points))
    assert total == certificate["total"]
    print(json.dumps({"status": "PASS", "points": total, "fresh_core_count": len(core),
                      "raw_product_maxima": maxima}, indent=2))


if __name__ == "__main__":
    main()
