# P48n body families and rational heads in dimensions 41–46

Checkpoint: **28 October 2026**. Six current constructions are supplied.

| Dimension | Total | Current certificate |
|---|---:|---|
| 41 | **1,460,036** | 1,460,034 body points and two pure-head poles |
| 42 | **1,991,445** | 1,991,439 body points and a pure-head hexagon |
| 43 | **2,919,782** | 118 complete families, 3,447 explicit tiny points and five pure heads |
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

Fresh integration replays passed **6,370,783** angle checks in 43D and
**2,427,690** in 45D. The portable verifiers reconstruct the parent embedding,
full rank-40 complement and gluing map, and check all 30,240 explicit shortest
bodies. Their complete-family populations retain the assumptions below.

## 43D continuous-head continuation

The current certificate has **2,919,782** points. Continuous head placement
after shared tiny-point deletions replaces 46 old points with 85, increasing
the explicit layer from 3,408 to **3,447** and the total by **39**.
The 118 complete families and five pure heads retain their previous inputs.

`partial43/verify.py` checks every angle using ambient 48-coordinate bodies.
`continuation43/verify_coordinates.py` independently reconstructs rank-40
coefficient coordinates and uses integer square-root cutoffs. Each checks
6,370,783 conditions, including 5,939,181 distinct tiny-point pairs and
22,001 short-difference refinements. `verify.py` runs both representations.

[The continuation notes](continuation43/README.md) give the continuous-polytope
method, exact coordinate definition and verification scope.
`continuation43/reproduce.py` rebuilds the final arrays from the retained
predecessor and rational delta, without a numerical search or optimizer.
`continuation43/search.py` and `continuation43/select.py` are optional further
search tools. `continuation43/negative_controls.py` checks that both verifiers
reject duplicated points and invalid pure-head angles.

The 24 October 2,919,743-point certificate is retained in
`continuation43/previous_points.npz`. Its 3,408-point layer had improved the
preceding complete tiny layer by 57. The older writeup target 2,911,990 and
its unavailable seven extension points remain historical; no current claim
depends on those missing files. See [the certification history](43_CERTIFICATION.md).

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
