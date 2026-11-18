"""Verify each pinned exchange and every pair in its final coordinate witness."""

import argparse
import importlib.util
import json
import os
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parent
REPOSITORY = ROOT.parents[2]
GEOMETRY = ROOT.parent / "dim28-takhanov-root-exchange" / "verify.py"
spec = importlib.util.spec_from_file_location("contact_sum_exchange_geometry", GEOMETRY)
geometry = importlib.util.module_from_spec(spec)
spec.loader.exec_module(geometry)


def integer_rows(points, count, dimension, label):
    geometry.need(points.shape == (count, dimension), f"{label} shape")
    geometry.need(points.dtype.kind in "iu", f"{label} integer dtype")
    geometry.need(np.all((points >= -11) & (points <= 11)), f"{label} coefficient bounds")
    return points.astype(np.int16)


def load_coordinates(path, count, dimension, key="points"):
    with np.load(path, allow_pickle=False) as archive:
        geometry.need(int(archive["squared_norm"]) == 128, "squared-norm metadata")
        geometry.need(int(archive["inner_product_limit"]) == 64, "inner-product metadata")
        geometry.need(int(archive["coordinate_denominator"]) == 2, "lattice scale")
        if key == "rational":
            geometry.need(
                archive["radical"].shape == (count, dimension)
                and not np.any(archive["radical"]),
                "source has zero radical part",
            )
        return integer_rows(archive[key], count, dimension, str(path))


def verify(record, direct=False, *, package=ROOT):
    dimension = record["dimension"]
    source = record["source"]
    source_path = REPOSITORY / source["path"]
    coordinate_path = package / record["coordinates"]
    exchange_path = package / record["exchange_file"]
    for path, digest in (
        (source_path, source["sha256"]),
        (coordinate_path, record["coordinates_sha256"]),
        (exchange_path, record["exchange_sha256"]),
    ):
        geometry.need(geometry.file_hash(path) == digest, f"file identity: {path.name}")

    incumbent = load_coordinates(source_path, source["count"], dimension, source["key"])
    final = load_coordinates(coordinate_path, record["count"], dimension)
    with np.load(exchange_path, allow_pickle=False) as archive:
        added = integer_rows(archive["added"], record["added"], dimension, "added rows")
        removed = integer_rows(archive["removed"], record["removed"], dimension, "removed rows")
    geometry.need(
        geometry.row_hash(final) == record["coordinate_array_sha256"],
        "canonical coordinate identity",
    )
    exchange = geometry.exchange_check(incumbent, final, removed, added)
    print("EXACT_EXCHANGE", dimension, json.dumps(exchange), flush=True)
    checked = geometry.integer_geometry(final, direct=direct)
    return {
        "verified": True,
        "dimension": dimension,
        "geometry": checked,
        "exchange": exchange,
        "coordinate_array_sha256": record["coordinate_array_sha256"],
    }


def main(package=ROOT):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--direct", action="store_true", help="use signed-int64 throughout")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    records = json.loads((package / "manifest.json").read_text())["records"]
    reports = []
    for record in records:
        report = verify(record, direct=args.direct, package=package)
        reports.append(report)
        print("VERIFIED", json.dumps(report), flush=True)
    if args.report is not None:
        args.report.write_text(json.dumps(reports, indent=2) + "\n")


if __name__ == "__main__":
    main()
