# Rationally deformed sign orbits in dimension 39

Checkpoint: **26 October 2026**. Exact lower bound: **769,180**.
Previous method: [Qihong Lin's dimension 39 construction](https://github.com/QH-Lin/kissing-numbers-33-39/tree/0c4e2bd486061b3fd55e128ac3497bfb12b73f07), **768,916** points.
The retained local root replacement construction has 768,924 points.

441,244 sparse points + 1,281 × 256 dense sign-orbit points = **769,180**.
Dense representatives are integer vectors normalized individually. Their
norms lie between 99,996 and 100,004. Complete Walsh enumeration checks
210,205,695 relative orbit products; the largest integer product is
4,997,124,511, below the allowed 4,999,600,008. Signed support bounds check
every dense point against every sparse family. All arithmetic is bounded
integer arithmetic with no tolerance. `source/proof.txt` supplies the proof.

```bash
python final.py
python verify.py
```

Python and NumPy suffice. The wrapper checks the complete SHA-256 manifest,
replays the full sparse-layer proof and every relative dense orbit pair in
a temporary directory, and checks the newly produced receipt. Parent data,
third-party attribution, both standalone verifiers, and original receipts
are included. No optimizer, search cache or download is required.

`source/verify_base.py` independently replays the inherited sparse layer.
