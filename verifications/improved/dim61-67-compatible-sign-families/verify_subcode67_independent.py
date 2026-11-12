"""Independent exhaustive sign audit for the 67D affine-subcode extension.

No discovery or production-verifier module is imported.  Binary kernels are
obtained by filtering all 65,536 possible words.  Every pair whose supports
intersect in more than eight positions is checked by direct integer products
of all its affine sign words, without the production parity-repair lemma.
Unchanged parts of the pinned published 67D construction remain inherited.
"""

from __future__ import annotations

import argparse
from collections import Counter
import copy
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np


ROOT = next(p for p in Path(__file__).resolve().parents
            if (p / "AGENTS.md").is_file() and (p / "verifications").is_dir())
SOURCE = ROOT / "verifications/improved/dim61-68-95-96-code-layers/source"
DEFAULT = Path(__file__).resolve().parent / "data/subcode67.json"
PREFIX = "verifications/improved/dim61-68-95-96-code-layers/source/"
PINS = {
    PREFIX + "shared_parity/expanded67_68.npz": "288ab6a9574e6163cabba45a0472204ec6a788300649a330edfc8ed9bce3140f",
    PREFIX + "dim61_96/code64_supports.npz": "cdbd5e5a7bbda9ddb84ee86052e2d4d00ccc331f0188de958a69ed83b0c529f6",
    PREFIX + "dim61_96/cyclic73.weight16.bin": "5d240938e26a7a477de62bc255439e743cb33120bb13b986a70c6b1fea0b5239",
    PREFIX + "dim61_96/shortening_support_search.json": "5422ed8f2b3be75267a57ea6fccc88f403a80e5bbe60dc64adc8edda51b61d88",
    PREFIX + "continuation49_72/netto.json": "95e62ec390f8ec100e8ddb9eb0bcb687d44fac9412457abcda00807e371a7b51",
}
ALL_WORDS = np.arange(65536, dtype=np.uint16)
WEIGHTS = np.fromiter((i.bit_count() for i in range(65536)), dtype=np.uint8)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def physical_positions(support):
    return tuple(i for i in range(67) if support >> i & 1)


@lru_cache(maxsize=None)
def enumerate_kernel(checks):
    """Brute-force parity filtering, independent of any row-reduction routine."""
    live = np.ones(65536, dtype=bool)
    for h in checks:
        live &= (WEIGHTS[ALL_WORDS & h] & 1) == 0
    return ALL_WORDS[live]


def decode_class(row, full=False):
    s = row["support"]
    h = row["checks"]
    t = row["shift"]
    require(type(s) is int and 0 < s < (1 << 67) and s.bit_count() == 16,
            "invalid weight-sixteen support")
    require(type(t) is int and 0 <= t < 65536, "invalid affine phase")
    require(type(h) is list and h and all(type(x) is int and 0 <= x < 65536 for x in h),
            "invalid binary checks")
    linear = enumerate_kernel(tuple(h))
    require(len(linear) > 1 and int(linear[0]) == 0, "empty or trivial kernel")
    minimum = int(WEIGHTS[linear[1:]].min())
    require(minimum >= 4, "minimum sign distance is below four")
    require(int(linear[-1]) == 65535, "class is not antipodal")
    if full:
        require(len(linear) == 2048, "replacement class lost full size")
    else:
        require(type(row["dimension"]) is int and 1 <= row["dimension"] <= 16,
                "invalid sign-code dimension")
        require(len(linear) == row["count"] == 1 << row["dimension"],
                "incorrect new class population")
    words = linear ^ t
    require(len(np.unique(words)) == len(words), "repeated affine sign word")
    return {"support": s, "positions": physical_positions(s), "words": words,
            "minimum_distance": minimum, "count": len(words), "checks": tuple(h),
            "shift": t}


def reconstruct_parent_supports():
    """Rebuild the 67-coordinate shortening directly from the pinned source."""
    data = np.fromfile(SOURCE / "dim61_96/cyclic73.weight16.bin", dtype="<u8")
    require(len(data) == 2 * 174324, "unexpected 73D support file length")
    choice = json.loads((SOURCE / "dim61_96/shortening_support_search.json").read_text())
    deleted = tuple(choice["dimensions"]["67"]["deleted_coordinates"])
    dual9 = (21, 37, 38, 49, 57, 63, 67, 70, 72)
    coords = tuple(i for i in range(1, 73) if i not in dual9) + (0,)
    coords += tuple(i for i in dual9 if i not in deleted)
    require(len(deleted) == len(set(deleted)) == 6 and set(deleted) <= set(dual9),
            "unexpected deletion coordinates")
    require(len(coords) == len(set(coords)) == 67 and not set(coords) & set(deleted),
            "invalid retained coordinate map")
    deletion_mask = sum(1 << i for i in deleted)
    supports = []
    for low, high in data.reshape(-1, 2):
        word = int(low) + (int(high) << 64)
        require(word.bit_length() <= 73 and word.bit_count() == 16, "invalid source support")
        if word & deletion_mask:
            continue
        support = sum(((word >> i) & 1) << j for j, i in enumerate(coords))
        supports.append(support)
    require(len(supports) == len(set(supports)) == 40428, "parent support count")
    with np.load(SOURCE / "dim61_96/code64_supports.npz", allow_pickle=False) as data64:
        embedded = [int(x) for x in data64["words"]]
    require(len(embedded) == len(set(embedded)) == 30828, "old 64D support count")
    require(set(embedded) == {s for s in supports if s < (1 << 64)}, "old embedding mismatch")
    return tuple(supports), deleted, coords


def source_sign_classes():
    old = {}
    new = {}
    with np.load(SOURCE / "shared_parity/expanded67_68.npz", allow_pickle=False) as data:
        for name, target in (("old", old), ("new", new)):
            supports = data[name + "_supports"]
            checks = data[name + "_duals"]
            shifts = data[name + "_shifts"]
            require(supports.shape == (1792,) and checks.shape == (1792, 5)
                    and shifts.shape == (1792,), "published sign family array shape")
            for s, h, t in zip(supports, checks, shifts):
                s = int(s)
                require(s not in target and 0 < s < (1 << 64) and s.bit_count() == 16,
                        "published support identity")
                target[s] = {"support": s, "checks": list(map(int, h)), "shift": int(t)}
    return old, new


def direct_pair(a, b):
    """Check every affine sign pair on common coordinates by exact products."""
    common = tuple(sorted(set(a["positions"]) & set(b["positions"])))
    ia = np.asarray([a["positions"].index(i) for i in common], dtype=np.uint16)
    ib = np.asarray([b["positions"].index(i) for i in common], dtype=np.uint16)
    # Values are +/-2; an int16 dot product has absolute sum <=64, far below
    # its overflow bound.  Disjoint support coordinates contribute exactly zero.
    x = 2 - 4 * ((a["words"][:, None] >> ia) & 1).astype(np.int16)
    y = 2 - 4 * ((b["words"][:, None] >> ib) & 1).astype(np.int16)
    products = x @ y.T
    maximum = int(products.max())
    if maximum > 32:
        i, j = np.unravel_index(int(products.argmax()), products.shape)
        raise ValueError("exact pair violation: " + json.dumps({
            "support_a": a["support"], "support_b": b["support"],
            "sign_word_a": int(a["words"][i]), "sign_word_b": int(b["words"][j]),
            "intersection": len(common), "inner_product": maximum}))
    return {"support_a": a["support"], "support_b": b["support"],
            "intersection": len(common), "sign_pairs": int(products.size),
            "maximum_inner_product": maximum}


def review67(certificate):
    raw = Path(certificate).read_bytes()
    c = json.loads(raw)
    require(c["format"] == "affine-sign-subcode-extension-v1" and c["dimension"] == 67,
            "certificate format")
    require(c["source_sha256"] == PINS and c["base_sign_family"] == "shared_parity/expanded67_68.npz",
            "certificate does not pin the expected published baseline")
    for name, expected in PINS.items():
        require(hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected,
                "source digest mismatch: " + name)
    parent, deleted, coords = reconstruct_parent_supports()
    parent_set = set(parent)
    fixed_old, fixed_new = source_sign_classes()
    require(set(fixed_old) <= parent_set and not set(fixed_new) & parent_set,
            "published modified and added support identities")

    changed = {}
    for row in c["modified_parent_classes"]:
        item = decode_class(row, full=True)
        s = item["support"]
        require(s not in changed and s in parent_set and s not in fixed_old,
                "replacement support is duplicate, absent, or previously modified")
        changed[s] = item
    added = [decode_class(row) for row in c["new_classes"]]
    require(added and len({a["support"] for a in added}) == len(added),
            "duplicate or absent new support")
    require(all(a["support"] < (1 << 64) and a["support"] not in parent_set
                and a["support"] not in fixed_new for a in added), "new support collision")

    # Existing full classes are replaced by equally large distance-four codes.
    # All their interactions with retained support classes are sign-independent.
    old_histogram = Counter()
    changed_old_pairs = 0
    for s in changed:
        for t in parent:
            if s == t:
                continue
            overlap = (s & t).bit_count()
            require(overlap <= 8, "changed old class meets another parent support above eight")
            old_histogram[overlap] += 1
            changed_old_pairs += 1
        for t in fixed_new:
            overlap = (s & t).bit_count()
            require(overlap <= 8, "changed old class meets a published added support above eight")
            old_histogram[overlap] += 1
            changed_old_pairs += 1

    decoded_fixed = {}
    new_histogram = Counter()
    exceptional = []
    low_overlap_sign_pairs = 0
    new_support_pairs = 0
    for index, a in enumerate(added):
        for label, supports in (("parent", parent), ("published_new", fixed_new)):
            for s in supports:
                overlap = (a["support"] & s).bit_count()
                new_histogram[overlap] += 1
                new_support_pairs += 1
                if overlap <= 8:
                    low_overlap_sign_pairs += a["count"] * 2048
                    continue
                if s in changed:
                    b = changed[s]
                else:
                    rows = fixed_old if label == "parent" else fixed_new
                    require(s in rows, "high-overlap pair has no specified retained sign code")
                    if s not in decoded_fixed:
                        decoded_fixed[s] = decode_class(rows[s], full=True)
                    b = decoded_fixed[s]
                record = direct_pair(a, b)
                record["kind"] = label
                exceptional.append(record)
        for b in added[:index]:
            overlap = (a["support"] & b["support"]).bit_count()
            new_histogram[overlap] += 1
            new_support_pairs += 1
            if overlap <= 8:
                low_overlap_sign_pairs += a["count"] * b["count"]
            else:
                record = direct_pair(a, b)
                record["kind"] = "new_new"
                exceptional.append(record)

    terms = json.loads((SOURCE / "continuation49_72/netto.json").read_text())["dimensions"]["67"]["terms"]
    expected_terms = {"dense64": 1 << 28, "weight16": 40428 * 2048,
                      "weight4": 11792 * 16, "poles": 134}
    require(terms == expected_terms, "unexpected inherited layer counts")
    base = sum(expected_terms.values()) + 1792 * 2048
    gain = sum(a["count"] for a in added)
    require(base == c["base_count"] == 355090822 and base + gain == c["total"],
            "incorrect total population")
    require(changed_old_pairs == len(changed) * (40428 - 1 + 1792), "changed pair accounting")
    require(new_support_pairs == len(added) * (40428 + 1792) + len(added) * (len(added) - 1) // 2,
            "new support pair accounting")
    # All sparse vectors have norm64; all cross-products with other layers
    # are bounded by amplitude alone, including every new or reordered class.
    require(16 * 2**2 == 64 and 16 * 2 * 1 <= 32 and 4 * 2 * 4 <= 32
            and 2 * 8 <= 32, "mixed-layer amplitude bound")
    return {
        "status": "PASS", "dimension": 67, "base_count": base, "gain": gain,
        "total": base + gain, "certificate_sha256": hashlib.sha256(raw).hexdigest(),
        "reviewer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "new_classes": len(added), "new_class_size_histogram": dict(Counter(a["count"] for a in added)),
        "new_minimum_distance_histogram": dict(Counter(a["minimum_distance"] for a in added)),
        "replaced_full_classes": len(changed), "replacement_population_each": 2048,
        "replacement_minimum_distance_histogram": dict(Counter(a["minimum_distance"] for a in changed.values())),
        "published_classes_bruteforce_decoded_for_affected_pairs": len(decoded_fixed),
        "changed_old_support_pairs": changed_old_pairs,
        "changed_old_overlap_histogram": dict(sorted(old_histogram.items())),
        "new_support_pairs": new_support_pairs,
        "new_overlap_histogram": dict(sorted(new_histogram.items())),
        "low_overlap_new_sign_pairs_certified_by_support": low_overlap_sign_pairs,
        "directly_enumerated_high_overlap_pairs": len(exceptional),
        "directly_enumerated_sign_pairs": sum(x["sign_pairs"] for x in exceptional),
        "maximum_direct_inner_product": max((x["maximum_inner_product"] for x in exceptional), default=None),
        "high_overlap_pairs": exceptional,
        "source_pins": PINS, "deleted_coordinates": deleted, "retained_coordinates": coords,
        "scope": "All changed codes and affected interactions are exact. Unchanged published baseline and Netto quadruples remain inherited; no full ambient point enumeration is claimed.",
    }


def rejection_controls(certificate):
    """Run malformed and geometrically false certificates through this checker."""
    raw = Path(certificate).read_bytes()
    original = json.loads(raw)
    outcomes = []
    cases = []

    bad = copy.deepcopy(original)
    bad["new_classes"][0]["shift"] ^= 2
    cases.append(("new_phase_bit_flip_preserves_population_and_norm", bad, "exact pair violation"))
    phase_bad = bad

    bad = copy.deepcopy(original)
    repair = next(row for row in bad["modified_parent_classes"]
                  if row["support"] == 1245535671217722628)
    repair["shift"] ^= 2
    cases.append(("old_phase_bit_flip_preserves_population_and_norm", bad, "exact pair violation"))

    bad = copy.deepcopy(original)
    row = bad["new_classes"][0]
    old_count = row["count"]
    row.update(checks=[65535], dimension=15, count=32768)
    bad["total"] += 32768 - old_count
    cases.append(("sign_code_distance_two_with_consistent_count", bad, "minimum sign distance"))

    bad = copy.deepcopy(original)
    bad["modified_parent_classes"][0]["checks"].append(3)
    cases.append(("replacement_full_code_loses_half_its_population", bad, "replacement class lost full size"))

    bad = copy.deepcopy(original)
    bad["new_classes"].append(copy.deepcopy(bad["new_classes"][0]))
    bad["total"] += bad["new_classes"][0]["count"]
    cases.append(("repeated_support_with_adjusted_claimed_count", bad, "duplicate or absent new support"))

    bad = copy.deepcopy(original)
    bad["total"] += 1
    cases.append(("inflated_total", bad, "incorrect total population"))

    bad = copy.deepcopy(original)
    key = next(iter(bad["source_sha256"]))
    bad["source_sha256"][key] = "0" * 64
    cases.append(("changed_baseline_source_pin", bad, "does not pin the expected published baseline"))

    with tempfile.TemporaryDirectory(prefix="kissing67_independent_controls_") as tmp:
        folder = Path(tmp)
        for name, changed, expected in cases:
            path = folder / (name + ".json")
            path.write_text(json.dumps(changed))
            try:
                review67(path)
            except ValueError as error:
                message = str(error)
                require(expected in message, name + " rejected for an unexpected reason: " + message)
                if name.startswith(("new_phase", "old_phase")):
                    require('"inner_product": 40' in message,
                            "phase control did not recover the expected concrete collision")
                outcomes.append({"control": name, "status": "REJECTED", "reason": message})
            else:
                raise ValueError("invalid certificate was accepted: " + name)

        path = folder / "optimized_python_bad_phase.json"
        path.write_text(json.dumps(phase_bad))
        run = subprocess.run([sys.executable, "-O", str(Path(__file__).resolve()),
                              "--certificate", str(path)], text=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60)
        require(run.returncode != 0 and "exact pair violation" in run.stderr
                and '"inner_product": 40' in run.stderr,
                "optimized Python failed to reject the concrete angle violation")
        outcomes.append({"control": "optimized_python_rejects_geometric_collision",
                         "status": "REJECTED", "returncode": run.returncode,
                         "reason": run.stderr.strip().splitlines()[-1]})

    return {"status": "PASS", "controls": outcomes, "control_count": len(outcomes),
            "certificate_sha256": hashlib.sha256(raw).hexdigest(),
            "reviewer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, default=DEFAULT)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--controls", action="store_true",
                        help="run eight deliberate corruption checks instead of the positive proof")
    args = parser.parse_args()
    result = rejection_controls(args.certificate) if args.controls else review67(args.certificate)
    if args.report:
        args.report.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k not in (
        "high_overlap_pairs", "source_pins", "retained_coordinates")}, indent=2))
