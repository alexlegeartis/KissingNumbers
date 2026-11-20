"""Reconstruct and exactly verify the 26D section replacement."""
import argparse
import hashlib
import gzip
import itertools
import json
import math
import os
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import numpy as np

HERE = Path(__file__).resolve().parent
I = [0, 2, 4, 5, 6, 10, 11, 23]
J = [i for i in range(24) if i not in I]
K = I + [24, 25]
SECTION_SHA256 = "5c49a45f7076f63c3bf4468e86ad20b64aa6b29a36c6d76715ede6df4a388551"
UPSTREAM_SHA256 = "8b704bb49503e3a4f6ad6b567690a7a2a9abfb30600f453313eb0510e6d29a94"


def leech_shell():
    """Standard Golay model of the 196,560 Leech roots, norm squared 128."""
    shifts = (0, 2, 4, 5, 6, 10, 11)
    gen = np.zeros((12, 24), dtype=np.int16)
    for i in range(12):
        gen[i, [(i + j) % 23 for j in shifts] + [23]] = 1
    msgs = (np.arange(4096)[:, None] >> np.arange(12)) & 1
    code = (msgs @ gen) % 2
    assert np.unique(code.sum(1), return_counts=True)[1].tolist() == [1, 759, 2576, 759, 1]
    pairs = np.zeros((1104, 24), dtype=np.int16)
    for k, (i, j, s, t) in enumerate(
        (i, j, s, t)
        for i, j in itertools.combinations(range(24), 2)
        for s, t in itertools.product((-8, 8), repeat=2)
    ):
        pairs[k, i], pairs[k, j] = s, t
    signs = 1 - 2 * ((np.arange(256)[:, None] >> np.arange(8)) & 1)
    signs = signs[((signs < 0).sum(1) % 2) == 0]
    octads = np.zeros((759, 128, 24), dtype=np.int16)
    for k, word in enumerate(code[code.sum(1) == 8]):
        octads[k][:, np.flatnonzero(word)] = 4 * signs
    odd = np.repeat(np.where(code, 2, -2)[:, None, :], 24, axis=1)
    odd[:, np.arange(24), np.arange(24)] *= -3
    a = np.concatenate([pairs, octads.reshape(-1, 24), odd.reshape(-1, 24)])
    assert a.shape == (196560, 24) and (np.sum(a * a, axis=1) == 128).all()
    return a


def exact_integer_pairs(points):
    # Binary32 GEMM is exact: entries are integers and every partial sum of
    # coordinate products has absolute value <= 128 by Cauchy-Schwarz.
    assert np.all((points.astype(np.int64) ** 2).sum(1) == 128)
    f = points.astype(np.float32)
    for i in range(0, len(f), 512):
        mat = f[i:i + 512] @ f[:i + 512].T
        r = np.arange(len(mat))
        mat[r, i + r] = -128
        assert mat.max() <= 64, f"integer conflict at row block {i}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("takhanov_d27", type=Path, nargs="?", help="optional pinned upstream d27_certificate.npz")
    parser.add_argument("--full", action="store_true", help="recheck every exterior pair instead of using the standard Leech-shell bound")
    args = parser.parse_args()
    assert __debug__, "Do not run this verifier with Python -O"
    packed = (HERE / "section.i8.gz").read_bytes()
    assert hashlib.sha256(packed).hexdigest() == SECTION_SHA256
    section = np.frombuffer(gzip.decompress(packed), dtype=np.int8).astype(np.int64).reshape(-1, 26)
    assert section.shape == (9028, 26)
    if args.takhanov_d27:
        assert hashlib.sha256(args.takhanov_d27.read_bytes()).hexdigest() == UPSTREAM_SHA256
        with np.load(args.takhanov_d27, allow_pickle=False) as src:
            ext = src["added_integer"].astype(np.int64)
        assert np.array_equal(section, np.delete(ext[ext[:, 25] == 0], 25, axis=1))
    exact_integer_pairs(section)
    shell = leech_shell()
    sf = section[:, :24].astype(np.float32)
    retained = []
    # Recompute every Leech-vs-section exclusion. No omitted-root list is trusted.
    for block in np.array_split(shell, 128):
        allowed = (sf @ block.astype(np.float32).T).max(0) <= 64
        retained.append(block[allowed])
    retained = np.concatenate(retained)
    assert retained.shape == (192144, 24)
    baseline = np.concatenate([section, np.pad(retained, ((0, 0), (0, 2)))])
    exterior = baseline[np.any(baseline[:, J] != 0, axis=1)]
    assert exterior.shape == (200864, 26)
    assert np.unique(exterior, axis=0).shape == exterior.shape
    if args.full:
        exact_integer_pairs(exterior)

    data = np.loadtxt(HERE / "kernel.txt", dtype=np.int64, ndmin=2)
    assert data.shape[1] == 11
    rows = [tuple(map(int, r)) for r in data]
    for i, a in enumerate(rows):
        assert a[-1] > 0 and math.gcd(*a) == 1
        assert sum(v * v for v in a[:-1]) == 128 * a[-1] ** 2
        for b in rows[:i]:
            dot = sum(x * y for x, y in zip(a[:-1], b[:-1]))
            assert dot <= 64 * a[-1] * b[-1]

    proj = np.unique(exterior[:, K], axis=0)
    assert proj.shape == (2601, 10)
    assert (
        10 * max(abs(int(x)) for x in proj.flat) * max(abs(int(x)) for x in data[:, :-1].flat) < 2**63
        and 64 * max(int(x) for x in data[:, -1]) < 2**63
    )
    assert np.all(proj @ data[:, :-1].T <= 64 * data[:, -1])
    print(json.dumps({
        "status": "PASS",
        "dimension": 26,
        "fixed_exterior": len(exterior),
        "kernel": len(rows),
        "rational_noninteger_kernel_rows": int((data[:, -1] > 1).sum()),
        "cardinality": len(exterior) + len(rows),
        "full_exterior_pair_check": args.full,
    }, indent=2))


if __name__ == "__main__":
    main()
