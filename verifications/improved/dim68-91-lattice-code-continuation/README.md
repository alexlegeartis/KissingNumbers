# Mixed-shell lattice sections and ternary-code head compression

These are construction lower bounds for ordinary kissing numbers, not lattice kissing numbers or global optimality statements. The comparison table uses attributable external-source lower bounds, with dimension embeddings and orthogonal sums explicitly labelled derived; it does not compare against our own preceding constructions. These are not a complete current worldwide-priority audit. All earlier certificates remain available. The manuscript is unchanged.

| Dimension | External-source baseline | New lower bound | Increase | External attribution |
|---:|---:|---:|---:|---|
| 68 | 331,737,984 | **472,598,888** | 42.461494% | [Edel–Rains–Sloane 64D](https://arxiv.org/abs/math/0207291), embedded *(derived)* |
| 69 | 331,737,984 | **780,388,888** | 135.242549% | [Edel–Rains–Sloane 64D](https://arxiv.org/abs/math/0207291), embedded *(derived)* |
| 71 | 331,737,984 | **2,623,578,936** | 690.858769% | [Edel–Rains–Sloane 64D](https://arxiv.org/abs/math/0207291), embedded *(derived)* |
| 87 | 6,218,178,164 | **7,000,788,678** | 12.585849% | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [15D tail](../../../common/published.py) *(derived)* |
| 88 | 6,218,179,920 | **9,329,710,944** | 50.039257% | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [16D tail](../../../common/published.py) *(derived)* |
| 89 | 6,218,181,330 | **12,433,498,768** | 99.953943% | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [17D tail](../../../common/published.py) *(derived)* |
| 90 | 6,218,183,254 | **20,961,198,018** | 237.095212% | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [18D tail](../../../common/published.py) *(derived)* |
| 91 | 6,218,187,548 | **20,961,198,378** | 237.094985% | [Nebe Γ72](https://arxiv.org/abs/1008.2862v3) ⊕ [19D tail](../../../common/published.py) *(derived)* |

## Methods and sources

**71D: normalized four-slice union across the norm-8 and norm-10 shells of an extremal 72-dimensional lattice.** The four slice populations and all ten pair inequalities are checked using rational arithmetic. [Construction and design argument](mixed_shell71_20261117/PROOF.md). Published inputs are the existence/minimum theorem for Nebe's Gamma72, its shell design theorem and extremal theta coefficients. See [Nebe, An even unimodular 72-dimensional lattice of minimum 8](https://arxiv.org/abs/1008.2862v3) and the norm-10 distributions in [Kominers, Proposition 3](https://arxiv.org/html/2606.08734). This improves the retained [Gamma72 central-section proof](../dim70-71-gamma72-cross-sections/).

**87–89D: antipodal ternary-simplex fibers with exact D3, D4 and D5 head assignments.** An explicit generator of the classical extended ternary quadratic-residue [48,24,15] code is checked for factorization, rank and self-duality. Selected six-coordinate prefix fibers each contain exactly 3^18 words. Frozen finite head assignments, every possible body-distance residue and sparse cross-layer bounds are checked exactly. Independent reports: [87/89D](ternary_heads87_89_independent_20261117/AUDIT.md), [88D](ternary_head88_independent_20261117/AUDIT.md).

**90–91D: three-coordinate shortening and antipodal triangle embedding, with the final RS9/Steiner sparse union below.** The length-45 code contributes 2·3^21 dense points. The91D dense construction appends a zero coordinate. [Independent generator and dense geometry proof, including the retained older sparse version](ternary_qr48_independent_20261117/AUDIT.md).

The code minimum 15 is a published premise, not inferred merely from rank or a computed weight enumerator: [Gabriele Nebe, On extremal self-dual ternary codes of length 48](https://doi.org/10.1155/2012/154281), Introduction and Lemma 3.5 ([author manuscript](https://arxiv.org/html/1109.1681)). The parent code and triangle embedding are classical ingredients. The finite compressed-head witnesses and their compatible unions are the constructions documented here.

The displaced 87D record is the [Gamma72 signed-tail construction](../dim73-88-gamma72-cosets/). The displaced 88–91D records are the [Z8 Golay sections and nonlinear sign fibers](../dim88-96-z8-golay/). These links attribute the actual retained constructions rather than equating them with a general literature bound.

## Reproduce

From this directory, run `python verify.py` (Python 3.10 or newer and SymPy). It invokes the frozen independent finite checkers, preserving their input hashes, and checks all eight declared totals. No solver, giant point listing, full code enumeration or tolerance-based angle test is needed. `final.py` supplies the registry values. Source minimum and design premises remain explicit; independent finite checks do not replace those theorems.

## Exact D4 fiber continuation in68–69D

The original frozen inputs and independent hash checks are preserved unchanged. See `d4_recovery_20261117/independent_audit/AUDIT_SCAN68_69.md`. All22 selected fibers are directly replayed with signed128-bit integer norm checks, saturated affine coordinates and exact rational pair inequalities. Floating pruning may omit points, so populations are certified lower counts, not exhaustive shell counts. These are spherical configurations, not lattice kissing configurations. The preceding extended-moment construction remains in `../dim68-69-gamma72-extended-moments/`.

## Checker inventory

`verify.py` runs the package. `verify_lattice.py` rebuilds the Leech completion table and invokes the unchanged `d4_recovery_20261117/independent_audit/replay_scan.py`. `mixed_shell71_20261117/build_table.py` builds the finite parity table; `mixed_shell71_20261117/verify71.py` checks design moments and all ten normalized pair cases. `ternary_qr48_independent_20261117/replay_generator.py` checks the classical generator and ranks; `ternary_qr48_independent_20261117/verify_sparse.py` checks sparse counts and geometry. `ternary_head88_independent_20261117/verify.py` and `ternary_heads87_89_independent_20261117/verify.py` preserve the frozen head checks. The lattice replay also requires NumPy and a C++17 compiler.

## Final sparse-layer replacement in87–91D

The final totals replace the historical XOR-zero four-support augmentation by three mutually compatible layers: extended Reed–Solomon graph supports of size10 with72 explicit nonlinear binary signs; all scaled D_n roots; and size8 supports obtained from a Steiner quadruple system on10 groups, near-perfect matchings and an XOR-zero label code, with128 even signs each. Do not union the old and new sparse versions. The72-sign code is not antipodally closed; do not add its negatives.

The complete original production and independent packages are preserved. Any historical “pending” wording in a production report records its earlier stage; final independent PASS reports are `rs9_sparse_independent_20261117/AUDIT.md`, `rs9_nonlinear72_independent_20261117/AUDIT.md`, and `rs9_sqs10_independent_20261117/AUDIT.md`. Their unchanged `verify.py` scripts are invoked by the package runner. Production regeneration scripts are `rs9_sparse_20261117/verify.py`, `rs9_nonlinear72_20261117/verify.py`, and `rs9_mixed_sqs10_20261117/verify.py`.

The explicit72-sign data are attributed to Östergård–Baicheva–Kolev's binary(10,72,3) code, [electronic data](https://users.aalto.fi/~pat/opt/72.txt), [author description](https://users.aalto.fi/~pat/72.html). The independent auditor could not retrieve that source URL; all2,556 distances of the explicit data are checked, so validity does not depend on bibliographic lookup. Reed–Solomon polynomial root counting and the displayed finite Steiner triple-incidence certificate provide the support intersection proofs. No new parent code or worldwide priority is claimed.

## Public Gamma72 provenance and completion theorem

The full derivation of all norm8 completion weights and their dependencies is in [the standalone theorem](public_gamma72_provenance_20261117/THEOREM.md), not just in earlier audit notes. `public_gamma72_provenance_20261117/verify.py` verifies the exact source isometry, all196,560 distinct norm4 witnesses, parity fiber counts and the complete table hash. The runner checks that these are the same source matrices used by the68/69 enumerators. The source is [Nebe–Sloane's Gamma72 catalogue entry](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/neb72.html); all5,184 Gram entries match its downloaded matrix. The downloaded source hash and canonical matrix hash are pinned in `source.json`. Source minimum8 and the standard rootless rank24 theta theorem remain explicit premises.
