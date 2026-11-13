"""Count incumbent module populations and potential blockers for all Golay octads."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from common.golay import golay24


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("source", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    z = np.load(args.source)
    X = z["points"].astype(np.int16)
    if int(z["squared_norm"]) == 128:
        assert not np.any(X % 2)
        X //= 2
    assert X.shape[1] == 30
    assert np.all(np.sum(X.astype(np.int64)**2, axis=1) == 32)
    words = golay24()
    octads = words[words.sum(axis=1) == 8]
    assert len(octads) == 759
    q = np.sum(X[:, 24:].astype(np.int64)**2, axis=1)
    shapes, counts = np.unique(np.c_[X[:, :24]**2, q], axis=0, return_counts=True)
    Q = shapes[:, 24]
    A = shapes[:, :24].astype(np.float32)
    records = []
    for start in range(0, len(octads), 32):
        H = A @ octads[start:start+32].astype(np.float32).T
        for offset, h in enumerate(H.T):
            I = np.flatnonzero(octads[start+offset])
            module = h+Q == 32
            external = (((Q == 8) & (h > 8)) | ((Q == 16) & (h > 0))) & ~module
            records.append(dict(octad=I.tolist(), module=int(counts[module].sum()),
                                equator=int(counts[module & (Q == 0)].sum()),
                                midpoints=int(counts[module & (Q == 8)].sum()),
                                heads=int(counts[module & (Q == 16)].sum()),
                                external_potential=int(counts[external].sum())))
    records.sort(key=lambda z: z["external_potential"])
    report = dict(source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
                  source_count=len(X), records=records,
                  scope="Projection-based potential obstruction counts, not actual replacement costs or an exclusion.")
    args.output.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(records[:3], indent=2))


if __name__ == "__main__":
    main()
