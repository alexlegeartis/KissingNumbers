"""Independently verify the complete 33-dimensional kissing configuration."""

import argparse
import hashlib
from itertools import combinations, product
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DIMENSION = 33


def require(condition, message):
    if not condition:
        raise ValueError(message)


def manifest():
    return json.loads((HERE / "manifest.json").read_text())


def load_witness():
    record = manifest()
    path = HERE / "witness.npz"
    require(hashlib.sha256(path.read_bytes()).hexdigest() == record["witness_sha256"],
            "Witness file digest mismatch")
    with np.load(path, allow_pickle=False) as archive:
        return record, dict(archive)


def indices(mask):
    return [axis for axis in range(DIMENSION) if mask >> axis & 1]


def verify(all_pairs=False):
    record, witness = load_witness()
    supports = list(map(int, witness["supports"]))
    parity = witness["parity"]
    require(witness["sparse"].dtype.kind in "iu", "Sparse rows must be integers")
    sparse = witness["sparse"].astype(np.int64)
    require(len(supports) == len(set(supports)), "Repeated complete octad support")
    require(all(0 < support < 1 << DIMENSION and support.bit_count() == 8
                for support in supports), "Invalid complete octad support")
    require(parity.shape == (len(supports),) and np.isin(parity, [0, 1]).all(),
            "Invalid octad parities")
    require(all((a & b).bit_count() <= 4 for a, b in combinations(supports, 2)),
            "Complete octads meet in more than four coordinates")

    require(sparse.ndim == 2 and sparse.shape[1] == DIMENSION, "Invalid sparse coordinates")
    require(np.all(np.abs(sparse) <= 5), "Unexpected sparse coordinate magnitude")
    require(np.all(np.einsum("ij,ij->i", sparse, sparse) == 32), "Invalid sparse squared norm")
    require(np.all(np.abs(sparse).sum(axis=1) <= 16), "Sparse point exceeds dense cross bound")
    require(len(np.unique(sparse, axis=0)) == len(sparse), "Repeated sparse point")
    sparse_pair_max = -32
    for start in range(0, len(sparse), 192):
        dots = sparse[start:start + 192] @ sparse.T
        dots[np.arange(len(dots)), start + np.arange(len(dots))] = -32
        sparse_pair_max = max(sparse_pair_max, int(dots.max()))
    require(sparse_pair_max <= 16, "Conflicting sparse points")

    mixed_max = -32
    for support, required_parity in zip(supports, parity):
        rows = sparse[:, indices(support)]
        abs_rows = np.abs(rows)
        maximum = 2 * abs_rows.sum(axis=1)
        wrong_parity = ((rows < 0).sum(axis=1) % 2 != required_parity)
        maximum -= 4 * abs_rows.min(axis=1) * wrong_parity
        mixed_max = max(mixed_max, int(maximum.max()))
    require(mixed_max <= 16, "Conflict with a complete signed octad")

    generators = list(map(int, witness["dense_generator"]))
    shift = int(witness["dense_shift"])
    positions = witness["dense_positions"].tolist()
    require(len(generators) == 17 and all(0 < value < 1 << 32 for value in generators)
            and 0 <= shift < 1 << 32, "Invalid dense code data")
    require(len(positions) == len(set(positions)) == 32
            and all(0 <= axis < DIMENSION for axis in positions), "Invalid dense embedding")
    differences = [0]
    for generator in generators:
        differences.extend([value ^ generator for value in differences])
    require(len(differences) == len(set(differences)), "Dependent dense generators")
    dense_distance = min(value.bit_count() for value in differences[1:])
    require(dense_distance >= 8, "Dense code distance is below eight")

    count = len(differences) + 128 * len(supports) + len(sparse)
    require(record["dimension"] == DIMENSION and record["count"] == count,
            "Manifest count mismatch")
    receipt = {"dimension": DIMENSION, "count": count, "current_gain": count - record["parent_count"],
               "dense_minimum_distance": dense_distance, "full_octads": len(supports),
               "explicit_points": len(sparse), "explicit_pairs": len(sparse) * (len(sparse) - 1) // 2,
               "maximum_explicit_inner_product": sparse_pair_max,
               "maximum_full_octad_to_explicit_inner_product": mixed_max}
    if all_pairs:
        signs = np.array(list(product([-2, 2], repeat=8)), dtype=np.int16)
        sign_parity = (signs < 0).sum(axis=1) % 2
        middle = np.zeros((128 * len(supports), DIMENSION), dtype=np.int16)
        for i, (support, required_parity) in enumerate(zip(supports, parity)):
            middle[128*i:128*(i+1), indices(support)] = signs[sign_parity == required_parity]
        dense = np.zeros((len(differences), DIMENSION), dtype=np.int16)
        words = np.asarray(differences, dtype=np.uint64) ^ np.uint64(shift)
        dense[:, positions] = 1 - 2*((words[:, None] >> np.arange(32, dtype="u8")) & 1).astype("i2")
        complete = np.concatenate((dense, middle, sparse.astype(np.int16))).astype(np.float32)
        # All products and intermediate sums have magnitude at most 33 * 5**2.
        # Thus the binary32 block Gram entries are exact integers.
        largest = -32
        for start in range(0, count, 256):
            block = complete[start:start + 256] @ complete[start:].T
            local = np.arange(len(block))
            block[local, local] = -32
            largest = max(largest, int(block.max()))
        require(largest <= 16, "Conflict in expanded geometry")
        receipt["maximum_inner_product_all_pairs"] = largest
        receipt["unordered_pairs"] = count * (count - 1) // 2
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--all-pairs", action="store_true")
    args = parser.parse_args()
    print(json.dumps(verify(args.all_pairs), indent=2))
