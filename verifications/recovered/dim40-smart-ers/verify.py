"""Check all Smart support pairs and exact inequalities for the three layers.

Pass a downloaded file as a positional argument, or explicitly use --download.
The 589824-word dense layer is a published theorem input, not reconstructed here.
"""
import argparse
import hashlib
from itertools import combinations
import json
from pathlib import Path
from urllib.request import urlopen

from final import results

def main():
    cert = json.loads(Path(__file__).with_name("certificate.json").read_text())
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("source_file", nargs="?", type=Path)
    mode.add_argument("--download", action="store_true")
    args = parser.parse_args()
    if args.download:
        with urlopen(cert["support_source"], timeout=20) as response:
            raw = response.read()
    else:
        raw = args.source_file.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == cert["support_sha256"]
    header, *words = raw.decode("ascii").splitlines()
    assert header == "$BASE=16" and all(len(s) == 10 for s in words)
    masks = [int(s, 16) for s in words]
    assert len(masks) == len(set(masks)) == cert["supports"]
    assert all(0 <= x < (1 << 40) and x.bit_count() == 8 for x in masks)
    pairs = 0
    for x, y in combinations(masks, 2):
        assert (x & y).bit_count() <= 4
        pairs += 1
    assert pairs == cert["checked_support_pairs"]
    # Squared norms 40, 8, 2. Each tuple is (upper dot, norm^2, norm^2).
    assert all(4*d*d <= a*b for d, a, b in
               [(20, 40, 40), (4, 8, 8), (1, 2, 2),
                (8, 40, 8), (2, 40, 2), (2, 8, 2)])
    print(json.dumps({**results()[0], "support_pairs": pairs,
                      "dense_layer": "published A(40,10) lower bound"}))

if __name__ == "__main__":
    main()
