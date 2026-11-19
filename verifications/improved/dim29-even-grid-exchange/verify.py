"""Exact integer replay for the 29D projected-grid exchange."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parent
BASE = ROOT.parent / "dim29-takhanov-auxiliary-exchange" / "data" / "coordinates29_232154.npz"
META = ROOT / "manifest.json"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def replay(base=BASE):
    m = json.loads(META.read_text())
    path = ROOT / "trade.npz"
    assert sha(path) == m["trade_sha256"]
    assert sha(base) == m["baseline_sha256"]
    with np.load(base, allow_pickle=False) as baseline:
        assert all(baseline[name].item() == value for name, value in (
            ("squared_norm", 128), ("inner_product_limit", 64), ("coordinate_denominator", 2)))
        A = baseline["points"].astype(np.int64)
    assert A.shape == (m["baseline_count"], 29)
    assert np.all(np.einsum("ij,ij->i", A, A) == 128)
    with np.load(path, allow_pickle=False) as trade:
        rays, norms, ids, cols = (trade[key] for key in ("rays", "squared_norms", "drop", "coordinates"))
    assert all(x.dtype.kind in "iu" for x in (rays, norms, ids, cols))
    assert rays.shape == (m["added"], cols.size)
    assert np.max(np.abs(rays)) <= 11
    assert cols.ndim == ids.ndim == norms.ndim == 1
    assert len(set(cols)) == len(cols) and np.all((0 <= cols) & (cols < 29))
    assert len(set(ids)) == len(ids) == m["removed"] and np.all((0 <= ids) & (ids < len(A)))
    assert len(norms) == len(rays) and np.all(np.einsum("ij,ij->i", rays, rays) == norms)
    scales = [next(k for k in range(1, 12) if k*k*int(n) == 128) for n in norms]
    U = np.zeros((len(rays), 29), dtype=np.int64)
    U[:, cols] = rays * np.array(scales)[:, None]
    assert np.all(np.einsum("ij,ij->i", U, U) == 128)
    assert len(set(map(tuple, U))) == len(U)
    D = U @ A.T
    assert D.max() < 128
    assert np.array_equal(np.flatnonzero((D > 64).any(axis=0)), np.sort(ids))
    dd = U @ U.T
    assert (np.triu(dd, 1) <= 64).all()
    keep = np.ones(len(A), dtype=bool)
    keep[ids] = False
    F = np.concatenate((A[keep], U)).astype(np.int16)
    assert F.shape == (m["count"], m["dimension"])
    assert hashlib.sha256(F.astype("<i2", copy=False).tobytes()).hexdigest() == m["final_little_endian_int16_sha256"]
    return F, {"dimension": 29, "total": len(F), "added": len(U), "removed": len(ids),
               "donor_incumbent_pairs": int(D.size),
               "donor_donor_pairs": len(U)*(len(U)-1)//2}

def all_pairs(F):
    # Every integer product and partial sum is exactly representable in binary32.
    assert 29 * int(np.max(np.abs(F)))**2 < 2**24
    W = F.astype(np.float32)
    pairs = 0
    for i in range(0, len(W), 768):
        B = W[i:i+768]
        s = B @ W[i:].T
        np.fill_diagonal(s[:, :len(B)], 0)
        assert s.max() <= 64
        pairs += int(B.shape[0])*(len(W)-i-len(B)) + len(B)*(len(B)-1)//2
    assert pairs == len(F)*(len(F)-1)//2
    return pairs

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, default=BASE)
    parser.add_argument("--full", action="store_true", help="check all final unordered pairs from scratch")
    parser.add_argument("--output", type=Path, help="export reconstructed coordinates to an NPZ")
    args = parser.parse_args()
    F, r = replay(args.baseline)
    if args.full:
        r["all_final_pairs"] = all_pairs(F)
    if args.output:
        np.savez_compressed(args.output, points=F, squared_norm=128,
                            inner_product_limit=64, coordinate_denominator=2)
    print(json.dumps({"status": "PASS", **r}, sort_keys=True))
