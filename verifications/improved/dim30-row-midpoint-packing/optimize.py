"""Recompute and solve a finite packing instance against the 248,416-point source."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ortools.sat.python import cp_model

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]


def load_instance(candidate_file):
    record, = json.loads((HERE / "manifest.json").read_text())["records"]
    parent = REPO / record["source"]["path"]
    if hashlib.sha256(parent.read_bytes()).hexdigest() != record["source"]["sha256"]:
        raise ValueError("Source file differs from the pinned construction")
    with np.load(parent, allow_pickle=False) as archive:
        source = archive["points"]
    pool = record["candidate_pool"]
    frozen_pool = candidate_file.resolve() == (HERE / pool["path"]).resolve()
    if frozen_pool and hashlib.sha256(candidate_file.read_bytes()).hexdigest() != pool["sha256"]:
        raise ValueError("Search pool differs from the pinned archive")
    with np.load(candidate_file, allow_pickle=False) as archive:
        candidates = archive["candidates"]
    if frozen_pool and len(candidates) != pool["count"]:
        raise ValueError("Search pool size differs from the manifest")
    for rows in (source, candidates):
        if rows.ndim != 2 or rows.shape[1] != 30 or rows.dtype.kind not in "iu":
            raise ValueError("Expected integer rows in 30 dimensions")
        if np.any(rows < -11) or np.any(rows > 11):
            raise ValueError("Unexpected coordinate magnitude")
        wide = rows.astype(np.int64)
        if np.any(np.einsum("ij,ij->i", wide, wide) != 128):
            raise ValueError("Expected squared norm 128")
    candidates = np.unique(candidates.astype(np.int16), axis=0)
    # Give each row its unique 60-byte key in both sets before finding the difference.
    def keys(rows):
        return np.ascontiguousarray(rows, dtype="<i2").view("V60").ravel()
    return source, candidates[~np.isin(keys(candidates), keys(source))]


def solve(source, candidates, workers, seconds):
    model = cp_model.CpModel()
    insert = [model.new_bool_var(f"insert_{i}") for i in range(len(candidates))]
    delete = {}
    source_matrix = source.T.astype(np.float32)
    for start in range(0, len(candidates), 64):
        # Every integer partial sum is bounded by 30 * 11**2, exact in binary32.
        conflicts = candidates[start:start + 64].astype(np.float32) @ source_matrix > 64
        for offset, row in enumerate(conflicts):
            i = start + offset
            for index in np.flatnonzero(row):
                index = int(index)
                if index not in delete:
                    delete[index] = model.new_bool_var(f"delete_{index}")
                model.add_implication(insert[i], delete[index])
            neighbors = np.flatnonzero(candidates[i + 1:] @ candidates[i] > 64) + i + 1
            if len(neighbors):
                model.add_bool_and([~insert[j] for j in neighbors]).only_enforce_if(insert[i])
    gain = sum(insert) - sum(delete.values())
    # Under equal gain, prefer a smaller exchange for the independent integer replay.
    model.maximize((len(candidates) + 1) * gain - sum(insert))
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = workers
    solver.parameters.max_time_in_seconds = seconds
    status = solver.solve(model)
    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        raise RuntimeError(f"Solver returned {solver.status_name(status)}")
    take = [i for i, variable in enumerate(insert) if solver.boolean_value(variable)]
    drop = [i for i, variable in delete.items() if solver.boolean_value(variable)]
    final = np.concatenate((np.delete(source, drop, axis=0), candidates[take]))
    print(json.dumps(dict(status=solver.status_name(status), candidates=len(candidates),
                          added=len(take), removed=len(drop), count=len(final))))
    return final


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, default=HERE / "data" / "candidates.npz")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--seconds", type=float, default=300)
    parser.add_argument("--output", type=Path, help="optional reconstructed coordinate archive")
    args = parser.parse_args()
    final = solve(*load_instance(args.candidates), args.workers, args.seconds)
    if args.output:
        final = final[np.lexsort(final.T[::-1])]
        np.savez_compressed(args.output, points=final, squared_norm=128,
                            inner_product_limit=64, coordinate_denominator=2)


if __name__ == "__main__":
    main()
