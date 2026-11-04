# New 96D support-pattern screen

No new lower bound resulted. This experiment expands row-occupancy patterns beyond the existing predominantly paired supports, while preserving the current 26022-support family.

The sufficient compatibility rule is intersection at most12 between any two weight24 supports. A support on at most12 of the24 four-cell rows also meets every transversal in at most12 places. This inexpensive support test deliberately ignores potentially useful sign-dependent compatibility; its failure is not a spherical impossibility theorem.

- Exhaustively tested all C(24,6)=134596 supports consisting of six complete rows: none satisfies the support rule against the existing family. `full_row_witnesses.bin` stores a blocker index for every combination in lexicographic order. `replay.py` independently checks the complete domain, source hash and exact integer overlaps. Minimum witnessed overlap14.
- Sampled1000000 supports with six triple rows and six singleton rows, deterministic mt19937_64 seed20261104. None satisfies the rule. This is a sampled negative result, not an exhaustive exclusion of that family. Standard-library shuffle implementations can differ; recorded compiler/runtime environment matters.
- Total screen8.473 seconds, child peak9344 KiB, 1358338904 exact popcount checks. `resources.txt` and `results.txt` record the run.

Replay the exhaustive certificate from the repository root:

```sh
python research/discovery_campaign_20261104/supports/replay.py
python -O research/discovery_campaign_20261104/supports/replay.py
```

Regenerate the binary input as little-endian uint64 interleaved low/high words from `verifications/improved/dim88-96-z8-golay/data/supports.npz`, compile `search.cpp` with existing `g++ -O3 -std=c++17`, and run the executable with input and witness-output paths. The executable expects exactly26022 source supports; the Python replay binds the actual source and certificate by SHA256. Expected source and certificate SHA256 values are checked explicitly, as are source weight, distinctness and encoding. Hashes detect a change to the selected evidence; they are not authenticated provenance. No existing record files were changed.
