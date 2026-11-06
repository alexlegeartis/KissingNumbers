# A 2,919,799-point kissing configuration in dimension 43

The configuration improves the repository's preceding 2,919,785 construction by fourteen points. It contains 2,916,330 points from 118 complete body families, 3,464 explicitly listed tiny-body points, and five pure heads. The complete families use an additive E8 Weyl discriminant relabeling of the earlier core; their individual populations change and are recomputed, while their total stays fixed. Continuous feasible-head searches then supply rational insertions. The relabeling is not asserted to lift to an isometry of the body lattice.

## Replay

Install NumPy and run `python verify.py` in this directory. No network access is required. `replay.json` records the result. The supplied run checks the parent embedding, full rank-40 complement, E8 gluing map, 30,240 explicit tiny vectors and their class counts, family populations against the stated identities, and every required angle. Head coordinates are rational integers divided by their stored denominator. All accepted inequalities use integer arithmetic; floating-point matrix inverses only propose matrices whose defining identities are then checked exactly.

`configuration_43.npz` specifies the complete families and pure heads. `explicit_points_43.npz` indexes the complete `tiny_bodies.npz` pool and gives rational heads. The realized vectors are the concatenated body/head vectors normalized to unit length. Body squared norms are `r3/3`; tiny bodies have squared norm 10/3. The integer raw tiny vectors are divided by sqrt(864). The head conventions are those of `verify.py`, with physical heads equal to stored rational heads divided by sqrt(3).

## Proof scope

This is a new exact geometry and angle replay, not a fresh enumeration of all 2.9 million vectors. It inherits the published minimum-six and 11-design properties of P48n and the complete-family population identities explained in `population_note.txt`. The representative exact normal-projection moment certificates are pinned in `inherited_e8_uniform_mass_certificate.json`; the original derivation script is included for provenance, not as a self-contained additional command.

Completeness of the tiny pool, needed for the tightened short-difference angle bounds, follows conditionally from the same population identities and the exhaustive parent-shell accounting:

1092000 + 240*116235 + 2160*9180 + 6720*495 + 9*30240 + 240 = 52416000.

The replay checks the 30,240 distinct tiny vectors, their actual gluing classes and norms. The population note explains the modular-form relation producing nonuniform class counts. These inherited premises are exposed rather than presented as a new full-shell enumeration.

`source_manifest.json` pins inherited source files and their hashes; `SHA256SUMS` pins the delivered files. `core_transform.npz` records the discriminant relabeling for reproduction; the angle proof does not require that transformation to lift to a body-lattice isometry. `comparison.json` records the accessible comparison sources and limitations: the reviewed repository bound was 2,919,785, while a live Cohn leaderboard retrieval failed. No exhaustive worldwide-priority assertion is made.

A second implementation is included: install NumPy and SymPy, then run `python independent_verify.py`. It reuses the supplied geometry replay and independently reconstructs coefficients and labels, then checks all raw-coordinate angles. Its report is written to `independent_replay/report.json`.

The final gain replaces tiny point index 3356 of the preceding 3,463-point tiny layer with two mutually compatible rational heads on tiny body 4137. `final_trade.json` records this shared-deletion step. All final angles are rechecked from the delivered point array; the preceding layer is not required to replay the proof.
