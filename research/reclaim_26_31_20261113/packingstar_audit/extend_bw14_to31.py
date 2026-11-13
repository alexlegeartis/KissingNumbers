"""Extend the successful BW14 orientation through every omitted tail axis."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
spec = importlib.util.spec_from_file_location("bw15", REPO / "research/reclaim_26_31_20261112/modules/bw15_transplant.py")
bw15 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bw15)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("source", type=Path)
    ap.add_argument("--module", type=Path, default=HERE / "bw14_module.npz")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    S = np.load(args.source)
    X = S["points"].astype(np.int16)
    if int(S["squared_norm"]) == 128:
        assert not np.any(X % 2)
        X //= 2
    assert X.shape == (251254, 31)
    I = bw15.NORMAL
    J = [j for j in range(24) if j not in I]
    proj = X[:, np.r_[I, np.arange(24, 31)]]
    module = ~np.any(X[:, J], axis=1)
    q = np.sum(X[:, 24:].astype(np.int64)**2, axis=1)
    h = np.sum(X[:, I].astype(np.int64)**2, axis=1)
    relevant = ((q == 8) & (h > 8)) | ((q == 16) & (h > 0))
    ids = np.flatnonzero(relevant)
    external = relevant & ~module
    oldset = {row.tobytes() for row in proj}
    R4 = np.load(args.module)["transformation_numerator"].astype(np.int32)
    base = bw15.bw15().astype(np.int32)
    records = []
    for omitted in range(7):
        keep = [k for k in range(7) if k != omitted]
        for sign in (-1, 1):
            Y = np.zeros_like(base)
            Y[:, :8] = base[:, :8] @ R4[:8, :8].T
            Y[:, np.array(keep)+8] = base[:, 8:14] @ R4[8:, 8:].T
            Y[:, 8+omitted] = sign * 4 * base[:, 14]
            assert not np.any(Y % 4)
            Y = (Y // 4).astype(np.int16)
            new = np.array([r for r in Y if r.tobytes() not in oldset], np.int16)
            graph = csr_matrix(new.astype(np.float32) @ proj[ids].astype(np.float32).T > 16)
            matching = maximum_bipartite_matching(graph, perm_type="column")
            gain = len(new) - int(np.sum(matching >= 0))
            cost = int(np.any(Y.astype(np.float32) @ proj[external].astype(np.float32).T > 16, axis=0).sum())
            np.savez_compressed(args.out / f"module_omit{omitted}_sign{sign}.npz", heads=Y)
            rec = dict(omitted=omitted, sign=sign, external_cost=cost, full_module_gain=58-cost, union_gain=gain)
            if gain > 0:
                reverse = np.full(len(ids), -1, np.int64)
                for i, j in enumerate(matching):
                    if j >= 0:
                        reverse[j] = i
                left = set(map(int, np.flatnonzero(matching < 0)))
                right, queue = set(), list(left)
                for i in queue:
                    for j0 in graph.indices[graph.indptr[i]:graph.indptr[i+1]]:
                        j = int(j0)
                        if matching[i] == j or j in right:
                            continue
                        right.add(j)
                        if reverse[j] >= 0 and int(reverse[j]) not in left:
                            left.add(int(reverse[j]))
                            queue.append(int(reverse[j]))
                remove = ids[sorted(right)]
                add15 = new[sorted(left)]
                added = np.zeros((len(add15), 31), np.int16)
                added[:, I], added[:, 24:] = add15[:, :8], add15[:, 8:]
                keepx = np.ones(len(X), bool)
                keepx[remove] = False
                F = np.r_[X[keepx], added]
                assert len(F) == len(X)+gain and len(np.unique(F, axis=0)) == len(F)
                for i, row in enumerate(added):
                    dots = np.zeros(len(F), np.int32)
                    for j in np.flatnonzero(row):
                        dots += F[:, j].astype(np.int32) * int(row[j])
                    dots[len(X)-len(remove)+i] = -32
                    assert dots.max() <= 16
                np.savez_compressed(args.out/f"exchange31_gain{gain}_omit{omitted}_sign{sign}.npz", removed=X[remove], added=added)
                rec.update(added=len(added), removed=len(remove), exact_changed_pairs="PASS")
            records.append(rec)
            print(json.dumps(rec), flush=True)
    (args.out / "report.json").write_text(json.dumps(records, indent=2)+"\n")


if __name__ == "__main__":
    main()
