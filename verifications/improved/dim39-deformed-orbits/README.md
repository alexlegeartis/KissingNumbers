# A 769,852-point kissing configuration in dimension 39

**Exact lower bound: τ(39) ≥ 769,852.** This is **32** above the preceding published 769,820 certificate. The unchanged sparse layer has 441,244 points. The dense layer has 10,269 representatives under a 32-element sign group, giving 328,608 points.

All representatives were jointly deformed; this is not an insertion into unchanged old coordinates. Both the primary integer proof and a separately implemented direct-integer audit pass. The complete dense check covers **1,687,391,811** relative products. The primary proof checks **35,161,056** actual middle-support bounds; the independent proof explicitly evaluates all **8,859,271,680** permitted middle sign sums.

## Replay

Python 3.10+, NumPy, and a C++17 compiler are required. Run without Python `-O`:

    python final.py
    OPENBLAS_NUM_THREADS=1 python verify.py

The package verifier checks the manifest, replays retained historical proofs, and replays the new primary and independent proofs in temporary directories. It also runs four independent extrema-witness checks and eight geometric rejection controls for the new certificate.

See `replacement29/PROOF.md`, `replacement29/DISCOVERY.json`, and `replacement29/independent/README.md`. The untouched original independent audit archive is `replacement29/independent_original.tar.gz`; its SHA256 is e1f83f77997aef9e6613fb0d9db4a2a70de090039b0a249f4d5f258f306e5431.

Certificate SHA256: dee846d08a9516ca12a55fcf89bde154e529ba6302d7146487bed92d18a58300.

## Comparison and provenance

The preceding published [769,820 certificate](https://github.com/alexlegeartis/KissingNumbers/commit/20723895cc4300b7e61d0cea980ecb63a00e3516) is retained in `replacement28/`. Earlier `source/`, `continuation/`, `continuation64/`, and `continuation32/` proofs also remain unchanged. The new certificate is in `replacement29/`. Lin's [primary paper](https://arxiv.org/abs/2611.03682) gives 768,916. These comparisons do not establish exhaustive worldwide priority. The manuscript remains unchanged.

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
- `replacement29/exact_middle.py`
- `replacement29/independent/check_witnesses.py`
- `replacement29/independent/negative_controls.py`
- `replacement29/independent/verify_independent.py`
- `replacement29/verify.py`
- `replacement29/verify_base.py`
- `source/verify.py`
- `source/verify_base.py`
- `verify.py`
