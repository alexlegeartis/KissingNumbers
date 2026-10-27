# Self-contained final coordinate replays

Preserve the sibling directory layout in FINAL_FULL_REPLAY_MANIFEST.json. Python >=3.9 plus NumPy, Linux resource module; no creator-code imports, network, baseline certificate, or recorded PASS report is needed.

Run sequentially:

    python research/27d_independent_verification_20261026/verify30_238796_full_final.py
    python research/27d_independent_verification_20261026/verify31_246596_full_final.py

Each recomputes the Golay shell and checks its exact code facts, all owner/norm/distinctness checks, ALL insertion pairs, free versus whole shell/insertions/free/pure, and ordinary pure pairs. SHELL_PROOF.md supplies the explicit algebraic shell pair lemma and midpoint argument, avoiding a 196560-square scan. Each then execs the exact exceptional-tail verifier, producing final counts238808 and246638. No inherited baseline proof is used in these full variants. Runtime may be several minutes each; execute serially. The30 template permits512MiB address space,31 permits256MiB; actual baseline usage was below256MiB.

These full final variants were source-adapted from already executed full320/383 checkers and syntax-checked, but NOT locally executed on final325/385. The final candidates were independently validated by pinned FULL-baseline incremental replays, retained separately. CI should run these full variants. Fresh outputs use *_full_replay.json and do not overwrite original incremental reports. Output elapsed time/RSS/report hashes naturally change on rerun; input and source pins remain immutable. Do not run Python with -O because assertions are the checks.
