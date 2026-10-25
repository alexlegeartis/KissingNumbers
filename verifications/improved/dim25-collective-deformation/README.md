# Collective deformation in dimension 25

Checkpoint: **26 October 2026**. Exact lower bound: **197,584**.
Previous method: Takhanov–Yun, [certificate and paper](https://github.com/k-nic/kissing_d25-31/tree/b3768310ba566da5c8947546266c103fbcdce181), **197,583** points.
The retained local lens construction has 197,579 points.

Starting from the 196,560-point Leech shell, delete 1,991 distinct owners,
insert 2,913 vertices on 971 circles, and add 102 rational moved points:
196,560 − 1,991 + 2,913 + 102 = **197,584**.
The certificate contains all coordinates and rotations; its precursor hash
records provenance and is not a replay dependency.

Circle coordinates lie in Q(sqrt(2),sqrt(3)). Exact rational half-angle
rotations preserve their squared norm 32. Circle feet lie on deleted Leech
endpoint segments, so convexity bounds their products with retained shell
points by 16. Distinct circles within one class have centre product at most
4 and radius sqrt(12), bounding every product by 16. Every remaining pair
is screened by outward integer intervals with 75-digit root enclosures;
inconclusive pairs use exact signed comparisons in Q(sqrt(2),sqrt(3)).
The rational moved points are checked against every retained Leech point.
The standard minimal-shell kissing property is the only lattice theorem used.

```bash
python final.py
python verify.py
```

Python and NumPy suffice. The wrapper checks `data/SHA256SUMS.txt`, rebuilds
the Leech shell from the extended binary Golay code using the supplied
`source/leech.py` and `source/golay.py`, and replays the complete proof in a
temporary directory. No cached shell, optimizer, network access or angular
tolerance is used. `source/collective_record25_verification.json` is the
original receipt; a new receipt is produced on each replay.
