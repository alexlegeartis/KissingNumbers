# A 769,820-point kissing configuration in dimension 39

**Exact lower bound: τ(39) ≥ 769,820.** This is 384 above the preceding published 769,436 certificate and 128 above the independently verified intermediate 769,692 construction. The unchanged sparse layer has 441,244 points. The dense layer has 10,268 representatives under a 32-element sign group, giving 328,576 points.

All representatives were jointly deformed; this is not an insertion into unchanged old coordinates. Both the primary integer proof and a separately implemented direct-integer audit pass. The complete dense check covers **1,687,063,204** relative products. The primary proof checks **35,157,632** actual middle-support bounds; the independent proof explicitly evaluates all permitted middle sign sums.

## Replay

Python 3.10+, NumPy, and a C++17 compiler are required. Run without Python `-O`:

    python final.py
    OPENBLAS_NUM_THREADS=1 python verify.py

The package verifier checks the manifest, replays retained historical proofs, and replays the new primary and independent proofs in temporary directories. It also runs independent extrema-witness checks and four geometric rejection controls.

See `replacement28/PROOF.md`, `replacement28/DISCOVERY.json`, and `replacement28/independent/README.md`. The untouched original independent audit archive is `replacement28/independent_original.tar.gz`; its SHA256 is b23c5695cc8f4055be9539da2966ee6070a8d1950e7b87c06d84f95541418fd5.

Certificate SHA256: dbc293625abdbc7f4377ee94e77344967877864ea7717482e34479ce794085df.

## Comparison and provenance

The previous published 769,436 certificate is retained in `continuation32/`. Earlier `source/`, `continuation/`, and `continuation64/` proofs remain unchanged. The new certificate is in `replacement28/`. Lin's [primary paper](https://arxiv.org/abs/2611.03682) gives 768,916. These comparisons do not establish exhaustive worldwide priority. The manuscript remains unchanged.

## Script inventory

The complete replay above invokes the current and retained proof routes. Helper scripts and their exact package paths are listed for auditability:

- `continuation/negative_controls.py`
- `continuation/verify.py`
- `continuation/verify_base.py`
- `continuation/verify_direct.py`
- `continuation32/exact_middle.py`
- `continuation32/independent/audit_sparse39.py`
- `continuation32/independent/audit_sparse_witnesses39.py`
- `continuation32/independent/prepare_dense.py`
- `continuation32/independent/replay.py`
- `continuation32/negative_controls.py`
- `continuation32/verify.py`
- `continuation32/verify_base.py`
- `continuation64/independent/audit.py`
- `continuation64/independent/negative_controls.py`
- `continuation64/verify.py`
- `continuation64/verify_base.py`
- `final.py`
- `replacement28/exact_middle.py`
- `replacement28/independent/check_witnesses.py`
- `replacement28/independent/negative_controls.py`
- `replacement28/independent/verify_independent.py`
- `replacement28/verify.py`
- `replacement28/verify_base.py`
- `source/verify.py`
- `source/verify_base.py`
- `verify.py`
