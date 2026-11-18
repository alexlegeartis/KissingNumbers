"""Verify ray-preserving replacement bounds for the three specified Leech codes."""

import argparse
import hashlib
import importlib.util
from itertools import combinations
import json
from pathlib import Path

import numpy as np

from verify_symmetric_lifts import evaluate_section24

ROOT = Path(__file__).resolve().parent
SOURCE_25 = ROOT.parent.parent / "verifications/improved/dim25-critical-cluster/source"
OCTAD = [0, 2, 4, 5, 6, 10, 11, 23]
COMPLEMENT = [i for i in range(24) if i not in OCTAD]
PUBLISHED_26_SHA256 = "0f0d9447269f6581057a4944308594a53d96b1f1310c80aedda20d7ecfc230f5"
ALTERNATE_26_SHA256 = "87b1e2465bb9232a9dde9da4f3c7d5cf06104dec3c5c89408cf138c7371b38c9"


def load_checker25():
    spec = importlib.util.spec_from_file_location(
        "critical_cluster_geometry", SOURCE_25 / "verify_general.py"
    )
    checker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checker)
    return checker


def extreme_gain(dimension):
    # Per sign, at most 2d nonpolar representatives and one pole.
    # At most one group of four children can use both poles.
    return 2 * (4 * dimension + 2) + 1


def bound_from_25(checker, shell, leech_cap):
    source = json.loads((SOURCE_25 / "certificate25.json").read_bytes())
    rows, denominators, removed, _, circles, *_ = checker.reconstruct(source, shell)
    circle_count = len(circles) // 3
    if len(circles) % 3 or not np.array_equal(
        circles, np.repeat(np.arange(circle_count), 3)
    ):
        raise ValueError("Circle structure differs from the critical-cluster witness")
    for i in range(circle_count):
        sample = range(3 * i, 3 * i + 3)
        for a, b in combinations(sample, 2):
            dot = checker.field_dot(rows[a], rows[b])
            dot[0] -= 15 * denominators[a] * denominators[b]
            if checker.field_sign(dot) > 0:
                break
        else:
            raise ValueError(f"Circle {i} has no pair above the threshold")
    general = len(rows) - len(circles)
    count = len(shell) - len(removed) + len(rows)
    if (count, circle_count, general) != (197586, 926, 247):
        raise ValueError("Unexpected 25D parent counts")
    return {
        "parent": count,
        "circle_triples_with_a_pair_above_15_over_32": circle_count,
        "balanced_non_Leech_cap": 2 * circle_count + general,
        "extreme_gain_cap": extreme_gain(25),
        "upper": count + leech_cap + 2 * circle_count + general + extreme_gain(25),
    }


def bound_from_26(path, expected_hash, shell, leech_cap):
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected_hash:
        raise ValueError(f"Incorrect parent certificate: {path.name}")
    with np.load(path, allow_pickle=False) as z:
        points = z["integer_points"].astype(np.int16)
        rational_count = len(z["fraction_denominators"])
    tail2 = np.sum(points[:, 24:].astype(np.int32) ** 2, axis=1)
    equatorial = points[tail2 == 0, :24]
    if any(np.any(row % 2) or tuple(row // 2) not in shell for row in equatorial):
        raise ValueError("Unexpected point in the equatorial Leech layer")
    head2 = np.sum(points[:, OCTAD].astype(np.int32) ** 2, axis=1)
    body2 = np.sum(points[:, COMPLEMENT].astype(np.int32) ** 2, axis=1)
    lifted = points[(head2 == 32) & (body2 == 64) & (tail2 == 32)]
    _, indices, sizes = np.unique(
        lifted[:, COMPLEMENT], axis=0, return_inverse=True, return_counts=True
    )
    if np.any(sizes != 8):
        raise ValueError("Incorrect square-lift fibre population")
    residual = (
        lifted[np.argsort(indices)][:, OCTAD + [24, 25]]
        .reshape(-1, 8, 10)
        .astype(np.int32)
    )
    gram = residual @ residual.transpose(0, 2, 1)
    diagonals = np.arange(8)
    if not np.all(gram[:, diagonals, diagonals] == 64):
        raise ValueError("Incorrect residual norm")
    gram[:, diagonals, diagonals] = 0
    if not np.all((gram == 0) | (gram == -64)) or not np.all(
        np.sum(gram == -64, axis=2) == 1
    ):
        raise ValueError("Fibre residuals do not form four signed axes")
    count = len(points) + rational_count
    remainder = count - len(equatorial) - len(lifted)
    if count != 201192:
        raise ValueError("Unexpected 26D parent count")
    return {
        "parent": count,
        "eight_point_fibres": len(sizes),
        "other_points": remainder,
        "balanced_non_Leech_cap": 2 * len(sizes) + remainder,
        "extreme_gain_cap": extreme_gain(26),
        "upper": count + leech_cap + 2 * len(sizes) + remainder + extreme_gain(26),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--published26", type=Path)
    args = parser.parse_args()
    # 2 sqrt(3) - 3 < 15/32 since (3 + 15/32)^2 > 12.
    if 111**2 <= 12 * 32**2:
        raise ValueError("Invalid rational threshold")
    checker = load_checker25()
    vectors, _ = checker.independent_shell()
    shell = set(vectors)
    leech_cap = int(evaluate_section24())
    result = {
        "Leech_parent_cap": leech_cap,
        "critical_cluster_25_to_26": bound_from_25(checker, shell, leech_cap),
        "alternate_26_to_27": bound_from_26(
            ROOT / "data/coordinates26_201192_alternate.npz",
            ALTERNATE_26_SHA256,
            shell,
            leech_cap,
        ),
    }
    if args.published26 is not None:
        result["published_26_to_27"] = bound_from_26(
            args.published26, PUBLISHED_26_SHA256, shell, leech_cap
        )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
