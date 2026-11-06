# A 769,436-point kissing configuration in dimension 39

**Exact lower bound: τ(39) ≥ 769,436**, an increase of **64** over the preceding 769,372 certificate and **256** over the initial 769,180 baseline. The 441,244-point sparse layer is unchanged. The new dense layer has 10,256 representatives under a 32-element sign group, giving 328,192 points. All dense representatives were jointly deformed; the count increase does not mean that the previous dense coordinates are retained.

Every dense point is the exact normalization of a supplied integer direction. The primary proof checks all **1,683,122,416** relative products and **35,116,544** actual support/representative bounds. The maximum integer product is **4,999,992,742,814,284**, below **4,999,999,600,000,008**. No optimizer output, numerical tolerance, external download, unverified lattice minimum or shell-population premise is trusted by the proof.

## Replay

Python 3.9+, NumPy and g++ with OpenMP are required. Run without `-O`:

    OPENBLAS_NUM_THREADS=1 python continuation32/verify.py
    python continuation32/independent/replay.py
    python final.py
    python verify.py

`continuation32/PROOF.md` gives the full proof. `continuation32/exact_middle.py` computes the exact middle-layer maximum using the actual support list and the sign group's projected parity. This replaces the stronger all-support bound used in preceding continuations. `continuation32/verify.py` checks bounded integer Walsh products and calls the retained sparse proof in `continuation32/verify_base.py`.

Certificate SHA256: e9f2dd59bf5169bf915be8e7d622a7ff14481471e2e0a96e0f9d90b46a5c6ad7.

## Comparison and provenance

The preceding 769,372 result is pinned at [commit b4c1fe9](https://github.com/alexlegeartis/KissingNumbers/commit/b4c1fe9b6339e1b90e12c6bb216908702f8aba90). Lin's [primary paper](https://arxiv.org/abs/2611.03682) gives 768,916. These comparisons do not establish exhaustive worldwide priority. The manuscript remains unchanged.

## Retained proofs and verification files

The original `source/` retains the 769,180 proof byte-for-byte, and `continuation/` retains 769,308. `continuation64/` retains 769,372, with its previous package README added for history. `final.py` derives the current count. The package `verify.py` checks the complete manifest and replays the retained and current proofs in isolated temporary directories.

Historical scripts are `source/verify.py`, `source/verify_base.py`, `continuation/verify.py`, `continuation/verify_base.py`, `continuation/verify_direct.py`, `continuation/negative_controls.py`, `continuation64/verify.py`, `continuation64/verify_base.py`, `continuation64/independent/audit.py` and `continuation64/independent/negative_controls.py`. The current `continuation32/verify.py`, `continuation32/verify_base.py` and `continuation32/exact_middle.py` independently derive the new bound from explicit data.

## Independent audit and rejection controls

`continuation32/independent/replay.py` regenerates the binary input using `prepare_dense.py`, compiles `dense_direct.cpp`, checks all dense pairs with direct integer arithmetic, and cross-checks 8,192 scalar products. `audit_sparse39.py` independently enumerates all actual support/sign patterns, checking 8,848,056,320 middle signed sums plus every other pair family. `audit_sparse_witnesses39.py` resolves the reported maxima to actual original point pairs using Python integers. No producer geometry implementation is imported by these audit scripts.

`continuation32/negative_controls.py` tests duplicate orbits, zero norms, duplicate group masks and a deliberate middle-layer collision. The independent replay also rejects duplicate and zero-norm binary records. Together with retained controls, the package replay checks sixteen corrupted inputs. The complete original independent audit is preserved in `continuation32/independent/original_audit.tar.gz`; all portable replay inputs and derived tables are regenerated rather than trusted from that archive.
