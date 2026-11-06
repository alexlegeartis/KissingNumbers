# A 769,372-point kissing configuration in dimension 39

**Exact verified lower bound: τ(39) ≥ 769,372**, another **64** points beyond the preceding 769,308 certificate and 192 beyond 769,180. This is a complete explicit configuration, with no inherited shell population or unverified minimum-norm premise. No exhaustive worldwide-priority claim is made.

The unchanged sparse layer has 441,244 points. The new dense layer has 5,127 representatives under a 64-element coordinate-sign group, giving 328,128 points. Every point is the exact unit normalization of a supplied integer row after a group action.

## Replay and proof

Run Python 3.9+ with NumPy and g++, without `-O`:

    OPENBLAS_NUM_THREADS=1 python continuation64/verify.py
    python continuation64/independent/audit.py
    python verify.py

`continuation64/verify.py` checks every 841,315,065 relative product with bounded int64 Walsh arithmetic. The independently written `continuation64/independent/audit.py` compiles `relative.cpp` and directly checks those products with exact integer arithmetic, independently reconstructing all sparse layers, group closure and normalized-ray distinctness. The maximum product 49,999,015,430,245 is below 49,999,960,000,008. Both proofs bind to certificate SHA256 bfea4089709e11e7b810da37f751868552ed0c39d34256c410e386531f5bc401.

The full argument is in `continuation64/PROOF.md` and the independent report in `continuation64/independent/AUDIT.md`. Only the independent checker's default path and explicit no-`-O` guard were adjusted for packaging; its original file checksums are preserved in `original_SHA256SUMS.txt`. Its two product rejection controls are recorded in `negative_controls.json`.

## Discovery and provenance

A further subgroup split breaks the 128-element symmetry of the preceding certificate into 64-element orbits. Keeping half of another orbit and jointly relaxing all representatives produces the new code. Full pair constraints are rebuilt and positive safety margins enforced before integer rounding. Neither optimizer output nor a search cache is trusted by either proof.

The baseline 769,308 is pinned at [commit 329642c](https://github.com/alexlegeartis/KissingNumbers/commit/329642cd3fb22f6d4fd8a33c6153dcd00a10c4da). Lin's [primary paper](https://arxiv.org/abs/2611.03682) gives 768,916. Both sources were freshly checked in the independent audit. The manuscript is unchanged.

## Retained history and verification files

`source/` retains 769,180 byte-for-byte. `continuation/` retains 769,308 byte-for-byte, with its previous package README added for provenance. `final.py` derives the current count. Package `verify.py` checks the manifest and replays all three constructions in isolated directories, including both proof routes and eight rejection controls for the preceding continuation, then both proof routes for the new continuation.

The historical `source/verify.py` and `source/verify_base.py`, the 128-group `continuation/verify.py`, `continuation/verify_base.py`, `continuation/verify_direct.py` and `continuation/negative_controls.py`, and the 64-group `continuation64/verify.py`, `continuation64/verify_base.py` and `continuation64/independent/audit.py` are all included and replayed as described above.

`continuation64/independent/negative_controls.py` recompiles the independent product checker and reproduces both rejection controls during every package replay.
