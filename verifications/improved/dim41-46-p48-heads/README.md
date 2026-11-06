# P48n body families and rational heads in dimensions 41–46

Checkpoint: **28 October 2026**. Six current constructions are supplied.

| Dimension | Total | Current certificate |
|---|---:|---|
| 41 | **1,460,036** | 1,460,034 body points and two pure-head poles |
| 42 | **1,991,445** | 1,991,439 body points and a pure-head hexagon |
| 43 | **2,919,799** | 118 relabeled complete families, 3,464 explicit tiny points and five pure heads |
| 44 | **4,939,281** | 1,034 complete body families |
| 45 | **7,655,056** | 2,187 complete body families and 16 pure heads |
| 46 | **13,387,797** | 13,387,770 body points and 27 pure heads |

All integer angle checks and head-dimension certificates replay. Complete-family
populations inherit the documented K40 identities and supplied tables; they
have **not all been freshly enumerated**. Explicit tiny bodies and the three
new 45D medium classes are checked directly. No optimality claim accompanies
these constructions. The manuscript remains its earlier snapshot.

The published comparison values in 47D (23,766,960) and 48D (52,416,000)
remain comparisons. They are not new results of this package.

## Replay

From this directory, with Python 3.9+ and NumPy:

```bash
python final.py
python verify.py
```

From the repository root:

```bash
python run_all.py --only p48-heads
```

`final.py` recomputes the six registered totals from the candidate arrays,
counting the explicit 43D tiny layer once. `verify.py` verifies the complete
manifest, exact integer angles, population formulas and head dimensions. It
runs the portable 43D and 45D verifiers in temporary copies to preserve the
supplied reports. It needs no optimizer, compiler, network or python-flint.
Run without `-O`.

The28October integration replays passed **6,370,783** angle checks in43D and
**2,427,690** in 45D. The portable verifiers reconstruct the parent embedding,
full rank-40 complement and gluing map, and check all 30,240 explicit shortest
bodies. Their complete-family populations retain the assumptions below.

## 43D finite body/head exchanges

The preceding certificate has **2,919,785** points. Three successive finite
exchanges improve the preceding2,919,782 by **three points**, leaving3,450
explicit tiny points. The118 complete families and five pure heads retain
the same population identities and P48n minimum/design inputs.

The separate43D task discovered the exchange mechanism. The generic
`search.body_head_exchange` operator independently generated the same final
point set from the30240 bodies and305 supplied head directions, without a
supplied delta or target count. This is seeded discovery, with no held-out
or worldwide-priority claim. The engine output is a proposal until full replay.

`partial43/verify.py` checks all angles in ambient48-coordinate bodies.
`continuation43/verify_coordinates.py` independently reconstructs the rank-40
coefficient metric and uses integer square-root cutoffs. Both certify the
current3450-point tiny layer. Neither freshly enumerates the large shells.

[The continuation notes](continuation43/README.md) describe the exact combined
delta and scope. `continuation43/reproduce.py` rebuilds the current arrays
from the retained3408-point predecessor. The former3447-point arrays, delta,
reports and notes are preserved in `continuation43/history_2919782/`.
`continuation43/negative_controls.py` rejects duplicates and bad pure angles.

Historical optional tools `continuation43/search.py` and
`continuation43/select_candidates.py` describe the earlier28October pool.
Their role and retained artifacts are documented in the continuation notes.

The older historical missing2,911,990-point writeup remains unnecessary;
see [the certification history](43_CERTIFICATION.md).

## Coordinated 45D movement

`coordinated45/verify.py` checks **7,655,056** points:

```text
1,084,590 + 4,624,200 + 168*9,180 + 778*495 + 18,900 + 16.
```

This is **4,455** above the previous installed 7,650,601. Coordinated motion
of large, medium and small heads retains every family of the 7,653,571-point
checkpoint and adds medium classes 77, 202 and 239. `added_bodies.npz` supplies
all 495 bodies of each new class. Their residues, norms, distinctness and
within-family distances replay. The sixteen pure heads remain in the same
exact rank-five head space. All angle checks use arbitrary-size integers.

## Construction and mathematical inputs

The pinned parent is the extremal even unimodular lattice P48n, of minimum
squared norm six, with 52,416,000 minimal vectors forming an 11-design.
These published facts are mathematical inputs. The package does not reprove
the parent's minimum. K40 is the integral orthogonal complement of an embedded
sqrt(3) E8. A complete body family has discriminant class c, squared norm
r/3, population w and head z/(D sqrt(3)); each resulting point is normalized
individually. The stored head coordinates lie in an exact subspace of rank
k, giving dimension 40+k. Pure heads have body zero.

For explicit bodies, `population_replay/data/e8_complement_reduced.npz`
supplies basis b and Gram q = b b^T/96. A body is `v b/(3 sqrt(96))`, with
v in its specified residue class and `v q v^T = 3r`. The portable READMEs and
`coordinated45/coordinates.py` give the normalized coordinate conventions.

Write `N_i = r_i D^2 + z_i.z_i`. The gluing map gives the body bound in units
`b_ij = r_i + r_j - 18 + mu(c_i-c_j)/4`. Coincident classes and radii use
`2r_i`; a diagonal family uses `2r_i - 18` for distinct bodies. With
`P_ij = b_ij D^2 + 2 z_i.z_j`, each sufficient angle condition is
`P_ij <= 0` or `P_ij^2 <= N_i N_j`. Pure-head conditions are checked in the
same way. The portable 43D replay additionally uses explicit body products
and the certified short-difference refinement.

`POPULATIONS.txt` and the portable `population_note.txt` derive populations
495 at norm four, 9,180 at 14/3, `116235 - 40*t_c` at 16/3 and 1,084,590
at norm six for the appropriate nonzero classes. Here t_c is the number of
norm-10/3 bodies in the class. The zero-class norm-six count is separately
1,092,000. The default replay checks these formulas and the supplied tiny
table; it does not freshly enumerate every inherited large shell.

## Optional historical population certificates

With NumPy, python-flint and a C++17 compiler named `g++`:

```bash
python verify.py --replay-populations --jobs 4
python certify43.py --jobs 4
```

The first command replays the older standalone 41D/42D/44D certificate in a
temporary copy, including 25 core classes. Its totals are 1,460,034,
1,989,954 and 4,939,281. The second freshly enumerates all 13 large shells
of the **older 2,878,155-point 43D seed**. Its successful recorded run is
`certification/recovered43.json`. Neither command freshly enumerates every
large shell in the latest six constructions.

## Files and provenance

- `bounds.json`: current ledger; `final.py` and `verify.py`: read-only integration drivers.
- `continuation43/`: preceding tiny layer, exact delta, independent coefficient verifier,
  reproduction, optional search tools and pinned receipts.
- `partial43/`: complete portable current 43D result, including `verify.py`,
  `verify_families.py`, explicit tiny points, geometry, population note and pinned reports.
- `coordinated45/`: complete portable 45D result, including `verify.py`,
  `verify_inherited.py`, `coordinates.py`, explicit new bodies and pinned reports.
- `data/candidate_41.npz`, `candidate_42.npz`, `candidate_44.npz`,
  `candidate_46.npz`: unchanged current certificates. `candidate_45.npz` is the
  retained historical 7,650,601-point certificate, superseded in the current ledger.
- `data/e8_geometry.npz`, `population_table.npz`, `pure_directions.json`:
  inherited geometry, tiny counts and the original pure-head construction.
- `data/old43_seed.npz`, `certify43.py`, `certification/recovered43.json`:
  the historical recovered 43D certificate and its fresh population replay.
- `POPULATIONS.txt`, `SOURCE_PROGRESS.txt`, `PROVENANCE.json`: unchanged
  original population derivation, fifth-pass report and initial source mappings.
- `UPDATE_PROVENANCE.json`: source paths, archive hash and every byte-preserving
  copy installed at the 24 October checkpoint.
- `population_replay/`: original standalone package, including `run_verify.py`,
  `geometry.py`, `verify_masses.py`, `old_e8_first_shells.py`, `exact_count.cpp`,
  `requirements.txt`, `README.txt`, original data, candidate arrays and reports.
- `data/SHA256SUMS.txt`: every installed file except itself, including the portable bundles.

The delivered package needs no access to the source research directories or ZIPs.
Finite sampled search optima and unsuccessful numerical candidates are not
unrestricted kissing-number upper bounds.

## 6 November 43D core relabeling

The [new portable certificate](weyl43/) gives **2,919,799** points (+14): 2,916,330 complete-family points, 3,464 explicit tiny points and five pure heads. The additive E8 discriminant relabeling changes actual classes and their individual populations; the family total is preserved after recomputation. The transform is not assumed to lift to a body-lattice isometry. Fresh parent/complement/gluing replay and two exact angle implementations pass. Large-shell populations, P48n minimum and design remain the documented inherited inputs. No fresh full-shell enumeration or exhaustive worldwide-priority claim is made.

The added replay scripts are `weyl43/verify.py`, `weyl43/verify_families.py` and `weyl43/independent_verify.py`. `weyl43/inherited_verify_masses.py` preserves the upstream moment-identity derivation as provenance; it is not a standalone command in this portable bundle.
