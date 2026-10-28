# Coxeter–Todd phase refills in dimensions 36 and 37

Checkpoint: **28 October 2026**. Both new bounds pass two exact proof representations.

| Dimension | Previous repository bound | This construction | Increase |
|---|---:|---:|---:|
| 36 | 501,672 | **501,687** | **15** |
| 37 | 505,488 | **505,503** | **15** |

The fixed base has 403,704 points. Replacing 147 retained original phase orbits
by 165 new ones gives 32,661 orbits, hence **501,687** equatorial points.
The 37D certificate selects 636 three-point latitude body orbits over 216 cells.
Its 1,908 bodies occur at both heights +2 and −2, giving **505,503** points.
The 15 added latitude orbits contribute 90 points; the improved equator supplies
another 54 points relative to the retained 505,359-point witness. All points have squared norm 12 and all
distinct pair products are at most six in the stated Gram basis.

These improve explicit repository witnesses. The historical Cohn comparisons
are 484,760 and 500,008, as preserved in
[the dated table](../../../common/cohn_20261024.json). The live Cohn retrieval
failed during the cited source review, so no fresh worldwide priority or optimality
claim is made. [The source review](../../../common/source_review_20261028_coxeter_todd.json)
records these comparisons. The manuscript remains at its documented preceding checkpoint.

## Reproduction

From this directory, with Python 3.9+ and NumPy:

```bash
python final.py
python verify.py
```

From the repository root:

```bash
python run_all.py --repository-only --only coxeter-todd
```

`final.py` derives the counts from the certificate populations. `verify.py`
checks the complete SHA-256 manifest and replays the retained phase witness,
then both new checkers and exact delta reconstruction in a temporary copy.
Installed certificates and receipts are preserved.

- `completion/verify.py`: standard-library rational short-shell enumeration,
  quotient-distance proofs, and all required relative-phase integer comparisons.
- `completion/verify_coordinates.py`: a separate rational enumerator and
  expanded integer-coordinate proof. It checks 97,983 phase points, 1,908
  latitude bodies, 61,667,370 phase pair products and all 1,819,278 distinct
  latitude-body pairs. Different-label lattice distances settle the other pairs.
- `completion/reproduce.py`: reconstructs the certificate byte for byte from
  the retained witness and the explicit accepted trades in `delta.json`.
- `completion/proof.md`: the complete angle argument, certificate convention,
  search logic and relation to PackingStar's correction/refill strategy.
- `phase/exact_lattice.py`: retained exact lattice enumerator used by the symbolic route.
- `phase/verify_sections.py`: optional historical moment and section obstruction.
- `completion/search/prepare.py`: prepares exact equator or latitude candidate pools.
- `completion/search/`: optional pool construction and deterministic-step
  finite-graph search. Search output is never accepted without both checkers.

Neither proof uses floating-point angle tolerances, an optimizer, a recorded
shell population or files outside this repository. The two new checkers share
no enumeration or angle code. Their common input is the certificate and its
stated mathematical lattice-distance argument.

## Why the search improved

The fixed lattice base allows 729,000 product cells, containing 6,561,000 phase
orbits. The previous equator occupied 3,979 cells. Search a larger candidate
pool, count shared blockers, and coordinate deletions and insertions rather
than requiring each candidate to fit by itself. The new equator uses 3,975
cells. For 37D, allow individual three-point orbits in an auxiliary body cell
instead of insisting on all nine points. This admits 15 additional orbits over
nine new cells, while retaining every old latitude orbit.

This is an independent finite-graph implementation of the correction/refill
idea described in [PackingStar](https://github.com/CDM1619/PackingStar/tree/67be454b83af4287e1a2b61e5bd8c3596a8061d6).
It does not reproduce or train their RL system. Their reported larger 26–31D
coordinates were not independently replayed here. Separate shared root-star
trades now improve 33–35D. Broader complete-support and signed-region searches
in this pass found no additional gain; this does not establish a ceiling.

## Retained history

`phase/` remains the byte-preserving 24 October witness: **501,633** in 36D
and **505,359** in 37D, with its original proof and exact replay. Its optional
`vectors36.npz` and `vectors37.npz` describe that historical witness, not the
new completion. Its optional moment calculation rules out improving 33–35D
by taking sections of that particular construction.

`source/` preserves the earlier complete-cell construction, **501,147** in
36D and **501,149** with two poles in 37D. `UPDATE_PROVENANCE.json` records its
original integration. Every new file is under `completion/`; the old source
and phase certificates, reports and coordinate exports remain unchanged.

Rows of all coefficient vectors use the supplied Gram basis. Unit Cartesian
coordinates are `coefficients @ cholesky(gram) / sqrt(12)`. The exact proofs
operate directly in the Gram basis. `data/SHA256SUMS.txt` covers the entire
package except itself, including the optional search sources.

The further refill removes 79 and adds 84 three-point phase orbits relative to
the preceding main certificate, for fifteen additional equatorial points.
All 636 latitude orbits are retained. The original phase/ certificate remains
501,633 / 505,359; the compact delta reproduces the current certificate from it.

`completion/negative_controls.py` tests both exact proof routes on duplicate
orbits and false counts. `completion/search/base_choices.py` exposes 756
independent base-cell replacements, verifies all old/new and new/new base
distances, and writes a search seed after compatible orbit retention. This
new search freedom has not yet produced a larger registered bound. The current
certificate uses the preceding fixed base. `completion/search/prepare.py`
also supports `--preserve-lifts` when preparing an equatorial search.
