"""Rebuild the core and verify the explicit angular inequalities in 40D."""

import argparse
import json
import os
from pathlib import Path
import sys
import tempfile

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "dim40-affine-centre-refill"))
from verify_expanded import AFFINE, enumerate_core, load_module


def verify(compiler, workers):
    assert __debug__, "Run with Python assertions enabled"
    data = json.loads((HERE / "certificate.json").read_text())
    gram = np.asarray(json.loads((HERE / "gram.json").read_text()), dtype=np.int64)
    residue = np.asarray(data["residue"], dtype=np.int64)
    assert gram.shape == (40, 40) and np.max(np.abs(gram)) < 2**20
    assert np.array_equal(gram, gram.T) and not np.any(np.diag(gram) % 2)
    assert residue.shape == (40,) and ((0 <= residue) & (residue < 3)).all()
    assert np.any(residue) and not np.any((residue @ gram) % 3)
    parent = load_module(AFFINE / "verify.py")
    assert parent.positive_definite(gram) == 6561
    enclosure = HERE / "enclosure.npz"
    with tempfile.TemporaryDirectory(prefix="kn40_p48q_") as scratch:
        work = Path(scratch)
        minimum = parent.check_minimum(gram, compiler, work, enclosure_path=enclosure)
        print("Rank-40 lattice minimum is at least six.", flush=True)
        core = enumerate_core(gram, residue, compiler, work, enclosure_path=enclosure,
                              workers=workers)
    assert np.all(core % 3 == residue)

    with np.load(HERE / "points.npz", allow_pickle=False) as source:
        groups = {54: core, **{norm: source[f"norm{norm}"] for norm in (36, 30, 6)}}
    for norm, rows in groups.items():
        assert rows.dtype == np.int8 and rows.shape == (data["counts"][str(norm)], 40)
        assert np.all(np.abs(rows.astype(np.int16)) < 128)
        assert len(np.unique(np.ascontiguousarray(rows).view("V40"))) == len(rows)
        for start in range(0, len(rows), 8192):
            a = rows[start:start + 8192].astype(np.int64)
            assert np.all((a @ gram * a).sum(1) == norm)
    print(f"Enumerated {len(core):,} core points and checked all vector norms.", flush=True)

    # Two different core rows have a nonzero lattice difference after division
    # by three. Its squared norm is at least six, so the core product is <= 27.
    maxima = {"54/54": 27}
    pairs = 0
    norms = list(groups)
    for i, r in enumerate(norms):
        for j in range(max(1, i), len(norms)):
            s = norms[j]
            maximum, _ = parent.max_product(groups[r], groups[s], gram, i == j, False)
            assert maximum <= 0 or 4 * maximum**2 <= r * s
            n, m = len(groups[r]), len(groups[s])
            pairs += n * (n - 1) // 2 if i == j else n * m
            maxima[f"{r}/{s}"] = maximum
            print(f"Squared norms {r} and {s}, maximum raw product {maximum}.", flush=True)

    total = sum(map(len, groups.values()))
    assert data["dimension"] == 40 and total == data["total"]
    return {"status": "PASS", "dimension": 40, "points": total, "core_count": len(core),
            "minimum_enumeration": minimum, "explicit_unordered_pairs": pairs,
            "raw_product_maxima": maxima}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--compiler", default="g++")
    parser.add_argument("--workers", type=int, default=8, help="parallel affine-core enumeration jobs")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    print(json.dumps(verify(args.compiler, args.workers), indent=2))
