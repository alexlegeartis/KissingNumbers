# Dimensions 41–46: independently normalized families over P48n

Construction checkpoint of **23 October 2026**, integrated from the two
`kissing41_48_research_part*.zip` archives. The external progress filename says
October 22; its contents describe the fifth pass on October 23 and supersede
older totals and READMEs in those archives.

## Bounds and verification scope

| Dimension | Working total | Supplied construction | Current replay |
|---|---:|---|---|
| 41 | 1,460,036 | 1,460,034 body points and two pure-head poles | Integer angles, exact rank and inherited populations |
| 42 | 1,991,445 | 1,991,439 body points and a pure-head hexagon | Integer angles, exact rank and inherited populations |
| 43 | 2,911,990 | **Latest certificate missing**; recovered older seed has 2,878,155 points | Older seed freshly certified; latest total excluded from the results registry |
| 44 | 4,939,281 | 1,034 complete body families | Integer angles and exact rank; an older standalone population replay is also supplied |
| 45 | 7,650,601 | 7,650,585 retained body points and 16 pure-head points | Integer angles, exact rank and inherited populations |
| 46 | 13,387,797 | 13,387,770 retained body points and 27 pure-head points | Integer angles, exact rank and inherited populations |
| 47 | 23,766,960 | Published comparison value | No new improvement claimed |
| 48 | 52,416,000 | Published parent-shell value | No new improvement claimed |

The five supplied current constructions enter the repository registry. This
does not promote a table lookup to a new ambient shell enumeration. The 43D
working total remains visible in `bounds.json` and this README, with its gap
explicit; neither the older seed nor a search log certifies the larger total.
No exhaustive literature audit or optimality claim accompanies this checkpoint.
The manuscript is an earlier snapshot and has not been updated for this package.

## Replay

From this directory, with Python 3.8+ and NumPy:

```bash
python final.py
python verify.py
```

`final.py` recomputes the five registered totals from the candidate arrays.
`verify.py` verifies the package manifest, inherited population formulas, every
internal and cross-family sufficient angle inequality, pure-head angles and
exact rational head ranks. It also checks the recovered older 43D seed, using
Python integers because its denominator is 10^10. It writes no repository files
and needs neither an optimizer nor python-flint. Run without `-O`.

From the repository root, the same consistency replay is registered as:

```bash
python run_all.py --only p48-heads
```

For the original population certificate, install the optional dependencies in
`population_replay/requirements.txt` (NumPy and python-flint) and provide a
C++17 compiler named `g++`, then run:

```bash
python verify.py --replay-populations --jobs 4
```

This runs `population_replay/run_verify.py` in a temporary copy. It replays the
published-parent embedding and design identities and enumerates the 25 core
classes used by the older 41D/42D/44D certificates. Their totals are 1,460,034,
1,989,954 and 4,939,281. It does **not** freshly enumerate every body class in
the latest five configurations. The latest population formulas additionally
use the analytical derivation in `POPULATIONS.txt`. The original reports are
pinned evidence of an earlier run, not evidence that this optional command has
run on the reader's machine. The default verifier always reports
`fresh_body_enumeration: false`.

### Fresh certification of the recovered 43D seed

```bash
python certify43.py --jobs 4
```

This independently reconstructs the parent embedding and gluing, replays the
exact 11-design population identities, freshly enumerates all 13 large shells
used by the recovered seed, and checks all 6,441 family inequalities using
arbitrary-size integers. It proves **2,878,155** points in dimension 43 without
the inherited modular population lemma. The command needs NumPy, python-flint
and `g++`; all generated files go into a temporary directory. Use
`--output NEW_DIRECTORY` to retain its reports. A successful local run is
recorded in `certification/recovered43.json`.

The added root writeup `kissing41-46_writeup.tex` specifies a different
245-family baseline of 2,911,983 points and seven additional points. It names
`kissing41-46_exact_certificates.zip` and `kissing43_extension_2911990.zip` as
the coordinate sources; neither is among the two supplied research archives.
The writeup explicitly delegates its finite proof inputs to these files and
does not give their coordinates. See [43D certificate status](43_CERTIFICATION.md).

## Construction

The parent is the pinned extremal even unimodular lattice P48n, with minimum
squared norm six and 52,416,000 minimal vectors. Its minimal shell is an
11-design. These published parent facts are assumptions; the supplied replay
does not reprove the 48D minimum.

Let K be the rank-40 integral orthogonal complement of an embedded sqrt(3) E8.
A row of a candidate specifies a body coset c in K*/K, an integer radius index
r (`r3`), a population w, and an eight-coordinate integer head z with common
denominator D. Take every body x of squared norm r/3 in that family, attach
the head z/(D sqrt(3)), and normalize the entire point individually. The
eight-coordinate heads have exact rank k, so the construction lies in 40+k
dimensions. Pure directions a are normalized independently and have body zero.

For explicit body coordinates, `population_replay/data/e8_complement_reduced.npz`
supplies the integral 40 by 48 basis b and Gram q = b b^T/96. A body has the
form x = v b/(3 sqrt(96)), with v in Z^40 in its specified congruence class
and v q v^T = 3r. `population_replay/geometry.py` reconstructs the full complement,
class labels and gluing witnesses. The compact labels in `data/e8_geometry.npz`
are E8 coordinates modulo three; `old_to_new` converts the older labels. The
30,240 norm-10/3 bodies are supplied explicitly in
`population_replay/data/dual40_min.npz`.

### Integer angle certificate

Write N_i = r_i D^2 + z_i.z_i. A gluing representative for a difference class
gives a body inner-product bound M_ij; in stored units b_ij = 6 M_ij is

    r_i + r_j - 18 + mu(c_i-c_j)/4.

Here `mu` is the squared length of a doubled E8 gluing vector. When two
different families have the same body class and radius, use b_ij = 2r_i to
include coincident bodies. On the diagonal use b_ii = 2r_i - 18, since
distinct bodies of one family differ by a lattice vector of squared norm at
least six. With P_ij = b_ij D^2 + 2 z_i.z_j, every sufficient family condition is

    P_ij <= 0, or P_ij^2 <= N_i N_j.

The pure/family condition is

    2 a.z_i <= 0, or (2 a.z_i)^2 <= ||a||^2 N_i,

and the pure/pure condition is the analogous inequality with the two direction
lengths. The verifier checks fixed-width arithmetic ranges before squaring;
the older 43D seed uses arbitrary-size integers throughout.

### Populations and their inheritance

The modular identity in `POPULATIONS.txt` derives, for the appropriate nonzero
quadratic types, the populations 495 at squared norm 4, 9,180 at 14/3,
116,235 - 40 t_c at 16/3, and 1,084,590 at norm six. Here t_c counts norm-10/3
bodies in the class. The zero-class norm-six population is separately
1,092,000. These two central populations must not be interchanged.

`data/population_table.npz` supplies t_c. The default verifier checks family
weights against these formulas and the supplied discriminant geometry; it
does not independently establish the modular-form lemma or regenerate t_c.
The original standalone replay verifies explicit small bodies, design moment
identities and selected core populations. See its `README.txt` for the complete
proof and the published parent references.

## What changed across the passes

41D uses the second-pass 1,460,034-point body construction with two added poles.
42D changes all eighteen 9,180-point families and reorganizes the medium layer
from 84 to 87 families before adding six pure heads. It is a joint replacement,
not three medium families appended to the old configuration.

44D retains the earlier 4,939,281-point construction. In the fifth pass, 45D
and 46D retain all their body families and add 16 and 27 pure directions in
their existing head spaces. In 45D the directions have squared length 80 and
normalized off-diagonal products -3/5 and 1/5; the largest positive squared
cross product with the retained points is 8/45. In 46D the corresponding
values are length 48, products -1/2 and 1/4, and cross bound 1/6. Both cross
bounds are strictly below 1/4.

The previously withdrawn 44D/45D moment cross-sections and the recovered
12,309,600-point 46D cross-section remain historical packages. Their status
does not forbid these different, larger constructions.

The source research also explores released outer families, actual tiny-point
exchanges, K41/K42/K45 body models and a ceiling for one complete-family 48D
model. Finite solver optima are not global kissing-number upper bounds. The
46D refined tiny-point pool has not had its full joint exchange problem
optimized. The source archives retain these searches; this package contains
the construction inputs and replay needed to review the bounds. No search
pool or optimizer is needed for the default verifier.

## Files and provenance

| Path | Role |
|---|---|
| `bounds.json` | Authoritative installed ledger, including missing-certificate and comparison statuses |
| `final.py` | Read-only totals driver used by the audit and maintainer workflow |
| `verify.py` | Manifest, population-table, exact rank and integer angle replay; optional older population replay |
| `certify43.py` | Fresh population and exact angle certificate for the recovered 43D seed |
| `certification/recovered43.json` | Recorded fresh 43D replay, shell counts and input hashes |
| `43_CERTIFICATION.md` | Proof scope and missing inputs for the larger 43D claim |
| `data/candidate_41.npz`, `data/candidate_42.npz` | Exact copies of the preferred third-pass candidates |
| `data/candidate_44.npz` | Exact copy of the standalone 44D certificate |
| `data/candidate_45.npz`, `data/candidate_46.npz` | Exact copies of the canonical fifth-pass extensions |
| `data/e8_geometry.npz`, `data/population_table.npz` | Inherited class geometry and short-body populations |
| `data/old43_seed.npz` | Older recovered 2,878,155-point seed; not the missing latest certificate |
| `data/pure_directions.json` | Original 45D/46D direction coordinates and construction description |
| `data/SHA256SUMS.txt` | Every installed source, document and data file except the manifest itself |
| `POPULATIONS.txt` | Original modular-form population derivation, unchanged |
| `SOURCE_PROGRESS.txt` | Latest authoritative fifth-pass report, unchanged |
| `PROVENANCE.json` | Source archive hashes and the source-to-installed mapping for every copied file |
| `population_replay/README.txt` | Original standalone explanation and replay instructions |
| `population_replay/run_verify.py` | Original geometry, design, core-count, rank and angle driver |
| `population_replay/geometry.py` | Original exact embedding, complement and gluing reconstruction |
| `population_replay/verify_masses.py` | Original exact design moment identities |
| `population_replay/old_e8_first_shells.py` | Original helper used by the population replay |
| `population_replay/exact_count.cpp` | Original enclosing-ellipsoid integer enumerator |
| `population_replay/data/` | Original parent lattice, complement, small bodies and moment certificates |
| `population_replay/certificate_*.npz` | Original older 41D/42D/44D configurations |
| `population_replay/reports/` | Pinned original reports and enumerator inputs |
| `population_replay/requirements.txt` | Optional standalone replay dependencies |

All 819 entries in the source checkpoint's SHA-256 manifest matched before
integration. Copied files preserve their original bytes. Python bytecode,
optimizer caches and exploratory pools are not installed. The delivered
construction package does not need access to the original ZIPs at runtime.
