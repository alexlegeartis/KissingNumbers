# A 769,308-point kissing configuration in dimension 39

**Exact verified lower bound: τ(39) ≥ 769,308**, an increase of **128** over the preceding repository construction of 769,180. Both the integer Walsh verifier and an independent direct-product verifier passed. This package does not establish worldwide priority.

The construction has **441,244 unchanged sparse points** and **2,563 dense sign orbits of order 128**, giving 328,064 dense points. Every dense point is the exact unit normalization of a supplied integer row after a coordinate-sign action. These integers, not the floating discovery coordinates, define the configuration.

The sparse layer and its source attribution are retained from the pinned [preceding package](https://github.com/alexlegeartis/KissingNumbers/tree/16eb556d7b491d5d819d2d538de269d68a87dbb4/verifications/improved/dim39-deformed-orbits). Its complete verification is rerun. There is no unverified lattice minimum, shell-population premise, optimizer, external download or search cache in the proof.

## Replay

Run with Python 3.9+ and NumPy, without `-O`:

    OPENBLAS_NUM_THREADS=1 python continuation/verify.py
    OPENBLAS_NUM_THREADS=1 python continuation/verify_direct.py

The standalone checker reconstructs the sparse layer, proves the sign group has order 128, checks all sparse/dense bounds, and evaluates all **420,575,485** required relative-orbit products using exact bounded int64 Walsh arithmetic.

The maximum integer product is **49,999,014,167,161**, below the permitted **49,999,960,000,008**. Every row has norm between 9,999,996 and 10,000,004. `continuation/PROOF.md` supplies the complete angle and count argument. `continuation/orbit_verification.json` is the replay receipt, not an assumed input.

Certificate SHA-256:

    522a6c8eec8edf8abc97a0b9c3bdd82191456b0d56f5f1f50fd44c94f3be0b45

## Comparison and provenance

The preceding repository count is 769,180. Qihong Lin's primary [paper](https://arxiv.org/html/2611.03682v1), checked on 6 November 2026, gives 768,916 in Theorem 1 and Table 2. It is the already-used comparison source, not a stronger newly superseding result. The arXiv page displays 22 October 2026 despite its 2611 identifier; this metadata inconsistency is recorded without resolving it. No comprehensive claim about all other current or unpublished constructions is made.

## Discovery

The existing dense base is invariant under a 256-element sign group contained in a 15-dimensional binary kernel. Its 128 kernel quotient shifts transport the successful nonbinary extra orbit. Superposing two transported deformation patterns produces a useful, initially infeasible seed. Reducing the retained sign symmetry to the subgroup fixing coordinate zero leaves 128 group elements and allows the dense representatives to move independently in the two halves.

The search kept one half of the additional orbit and jointly optimized all 2,563 representatives. It rebuilt the full relative-pair constraint list after optimization. Safety margins were then imposed on both dense products and sparse-layer bounds before rounding to scale 10,000,000. The first rounding attempt failed a sparse bound and was discarded; only the refined certificate named above is accepted.

`continuation/DISCOVERY.json` records the construction chain. Discovery optimization is not part of verification and its numerical output is not trusted by either proof.

## Retained history and package replay

The unchanged `source/` directory retains the 769,180-point certificate and its original proof. `final.py` derives the new count. The package `verify.py` checks the complete manifest, replays the historical certificate, then both new proof routes and their rejection controls in isolated temporary directories. No file under `paper/` is changed.

    python final.py
    python verify.py

### Verification files

`source/verify_base.py` retains the historical sparse-layer replay. `continuation/verify_base.py` replays that same unchanged sparse input for the new certificate. `continuation/negative_controls.py` requires both new proof routes to reject duplicate orbits, false totals, zero-norm rows and duplicate group elements.
