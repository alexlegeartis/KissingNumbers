"""Freshly enumerate the affine core and verify the explicit cross-layer pairs."""

import argparse
from fractions import Fraction
import json
import os
from pathlib import Path
import subprocess
import tempfile

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")

import numpy as np

from verify import AFFINE, CORE_COUNT, GLUING, HERE, load_module, witness


def enumerate_core(gram, residue, compiler, work):
    enclosure = np.load(AFFINE / "enclosure.npz", allow_pickle=False)
    scale = int(enclosure["P"])
    diagonal = enclosure["D"].astype(object)
    lower = enclosure["M"].astype(object)
    inverse = enclosure["inverse_times_three"].astype(object)
    assert np.array_equal(gram.astype(object) @ inverse, 3 * np.eye(40, dtype=object))
    assert np.all(np.diag(lower) == scale) and not np.any(np.triu(lower, 1))
    assert np.all(diagonal > 0)
    error = scale**3 * gram.astype(object) - (lower * diagonal) @ lower.T
    eta = max(sum(abs(int(x)) for x in row) for row in error)
    trace = Fraction(sum(int(x) for x in np.diag(inverse)), 3)
    # If xQx=54, then ||x||^2 <= 54 trace(Q^-1). The error
    # between the rational LDL form and Q has operator norm <= eta.
    bound = 54 * (scale**3 + eta * trace)
    radius = -(-bound.numerator // bound.denominator)
    data = work / "core.input"
    with data.open("w", encoding="ascii") as stream:
        print(scale, radius, file=stream)
        for row in (diagonal, *lower, *gram, residue):
            print(*row, file=stream)
    binary = work / "enumerate_core"
    result = work / "core.i8"
    subprocess.run([compiler, "-O3", "-std=c++17", str(HERE / "enumerate_core.cpp"),
                    "-o", str(binary)], check=True)
    out = subprocess.check_output([str(binary), str(data), str(result)], text=True)
    core = np.fromfile(result, dtype=np.int8).reshape(-1, 40)
    assert int(out) == len(core) == CORE_COUNT
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
