# P48 continuation, independent code review, and a 67D affine subcode extension

Checkpoint: 13 November 2026. Discovery work is confined to this research
directory. Final witnesses and independent checkers are copied into the
canonical proof package by the campaign integration step.

## Outcome

The final affine sign-subcode witness gives **355,095,446 points in dimension 67**,
adding **4,624** to the retained 355,090,822-point construction. The new geometric
idea is to allow smaller affine sign codes on additional weight-sixteen
supports. They can obey more parity constraints than a full extended Hamming
code. The old classes affected by these constraints retain all 2,048 signs.

The final certificate is `subcode67_joint_certificate.json`. Its standalone
exact checker is `verify_subcode67.py`; `subcode67_joint_verification.json`
records the complete local replay. The first, separately preserved certificate
`subcode67_certificate.json` adds 512 points by the same mechanism.
The integrated publication package is
[`dim61-67-compatible-sign-families`](../../../verifications/improved/dim61-67-compatible-sign-families/README.md).

The dimensions 40–48 searches below found **no improvement**. Numerical search
outputs, finite-pool exhaustion and unsuccessful relaxations are not general
upper bounds on kissing numbers.

## The 67D construction

The unchanged parent has 40,428 ordinary weight-sixteen supports from its
73-coordinate code shortening, 1,792 further signed supports, dense vectors,
Netto quadruple signs, and axial points. All vectors have squared norm 64,
and the allowed inner product is 32.

A new support carries 16 entries of magnitude 2. Its signs form an affine binary
code `shift + kernel(checks)`. Direct enumeration verifies the cardinality and
minimum Hamming distance of every new code. The final 49 classes have sizes:

| Signs per class | Classes | New points |
|---:|---:|---:|
|512|1|512|
|256|6|1,536|
|128|7|896|
|64|19|1,216|
|32|13|416|
|16|3|48|
|**Total**|**49**|**4,624**|

When two supports meet in at most eight positions, their inner product is at
most 32 for every choice of signs. A nine- or ten-position intersection is
repaired by a common weight-eight parity check with opposite affine phases.
This forces a sign disagreement in the intersection, giving product at most
`4*(intersection_size-2) <= 32`.

The joint construction changes 357 previously untouched parent sign classes:
119 on the original 64 coordinates and 238 with nonzero tail coordinates.
For every shared old class, all demanded checks are jointly completed to an
RM(1,4) dual, and their affine phases are solved consistently. Thus each changed
old class still has 2,048 signs and minimum distance four. All new/new support
intersections are at most eight. The checker verifies every changed-old versus
parent interaction, every new versus parent or published-new interaction, and
all new/new interactions using integers.

The unchanged 67D base construction and its Netto quadruple design remain
explicit inputs. The integrated package separately replays the full original
1,792-block baseline and every retained parent-support pair, so a minimum-distance
claim for the entire 73-bit mother code is unnecessary for these support layers.
The new changes use no floating-point angle tolerance.

### Reproduce the finite proof

From the repository root, on an available reserved CPU allocation:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/dim32plus_continue_20261113/middim/verify_subcode67.py research/dim32plus_continue_20261113/middim/subcode67_joint_certificate.json
```

The checker reads only this certificate and pinned files from existing proof
packages. Discovery pools and compiler products are not proof inputs.
The independent audit separately enumerates all sign patterns and directly
checks 115,998,720 products across the 758 exceptional support pairs. Its
portable copy is `verify_subcode67_independent.py` in the canonical package.

### Discovery scripts

`prepare_subcode67.py` reconstructs the fixed 1,792-block family and full 67D
parent supports. `subcode67_search.cpp` scans the supplied 112,896 candidate
supports, solves forced affine phases, and intersects a Hamming code with the
necessary checks. It found 67 phase-feasible candidates. `freeze_subcode67.py`
exports the initial 512-point witness. `subcode67_collect.cpp` retains the best
subcode found for every phase-feasible support. `joint_subcode67.py` then performs
a bounded greedy search with strict new/new overlap at most eight, exact old-dual
completion, and consistent old phases. The final run lasted 60.00 seconds and
6,545 iterations; its best witness was found at iteration 3. These are finite
search outcomes, with no optimality claim.

## Independent review of the broader code-family results

`review_codes.py` imports no code from the campaign's `codes` verifiers. It
constructs sign kernels by column elimination, checks every changed support
pair using packed Python integers, explicitly enumerates affine signs for the
selected 63D section, and independently proves the 61D/62D compression inequalities
using rational enclosures of square roots. It also checks the 65D cap transfer.

The initial review is `independent_code_review.json`. Passing a certificate path
repeats the review for that exact file and writes a separately named receipt.
For a global family it selects and independently recounts a valid 63D octad
section; it does not require a supplied optimizer score or section receipt.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/dim32plus_continue_20261113/middim/review_codes.py research/dim32plus_continue_20261113/codes/global_final_certificate.json
```

The review found no mathematical gap. It identified two choices that must be
stated explicitly: untouched old Hamming codes may be reordered independently
for the selected 63D octad, and the 65D owner construction uses the inherited BCH
dense code in place of the 64D dense code. Both are valid because the relevant
cross-layer bounds depend only on amplitudes or support intersections. Standalone
assertion-based entry points were also identified for explicit `-O` rejection.
The root package implements those guards. The portable independent reviewer
is `verify_independent.py` in the canonical package; it imports no research
code, exposes a read-only `check()` function, and has a deterministic receipt
in `independent_verification.json`. The final global family has 1,940 new
supports and 303 repaired new/new pairs, and its independently recounted
63D section contains 184,393,988 points.

## Finite negative searches in 40–48D

### Whole-profile compression

`adaptive_profiles.py` and `twohead_profiles.py` aggregate the actual P48n
normal fibers in the existing sqrt(3)E8 complement package. They replace a
normal projection by one or two adjustable head coordinates, with a separate
body scale for each profile. The following are **numerical discovery scores**;
none is promoted to a newly certified kissing bound.

| Normal subsystem and output heads | Dimension | Best finite score | Bound to beat |
|---|---:|---:|---:|
|A2, one head|47|23,766,960|23,766,960|
|A3, two heads|47|22,216,330|23,766,960|
|A3, one head|46|12,309,600|13,387,797|
|D4, two heads|46|12,308,040|13,387,797|
|D4, one head|45|6,984,558|7,655,056|
|D5, one head|44|4,149,492|4,939,281|
|E6, one head|43|2,550,816|2,919,803|
|E7, one head|42|1,791,392|1,991,445|

The larger subsystem runs omit profiles of population below 1,000. No search
is exhaustive. `head_search.cpp` and `multihead_search.cpp` are compatibility-
graph discovery engines, not proof verifiers.

`normal_models.py` explores two different normal lattices. An actual P48p
three-vector embedding with Gram diagonal 6 and off-diagonal -2 is retained in
`bcc3_embedding.json`; its normal lattice is a scaled body-centred cubic lattice.
The one-head model scored 11,306,800, below the 46D record. The orthogonal-frame
model's moment relaxation was too loose to provide a competitive lower bound.
`pair_normals.py` explores overlaps -1 and -2 and recovers the known 47D section.

A complete-family 40D search, `family40_search.py` and `.cpp`, deletes the entire
norm-six core and searches smaller-radius body families. Its best finite score
was 929,940. Exact universal family-angle tests pass for that finite selection,
but it is below the current 1,126,314 and is not a new result. An additional
10,000-step E8 Weyl relabeling test changed none of the complete-family totals
in 41D, 42D, 44D, 45D, or 46D; `weyl_reweight.json` records those observed ranges.

### Joint medium/tiny exchanges in 43D

`medium43_pool.py` builds exact rational proposals and exact blocking-point sets
from the current medium-splitting certificate. The initial 5,898 proposals yielded
978 two-copy split options. A 120-second mixed-integer search found gain zero
with remaining finite-model upper bound five; this is an unresolved finite model,
not an impossibility proof.

`medium43_head_vertices.py` generated 81 further rational head templates from
the fixed-major-family polytope. The enlarged pool has 25,515 proposals.
`medium43_small_exchange.py` found no one-for-two or two-for-three exchange.
`prepare_medium43_options.py` reconstructs the option tables; the 44 low-deletion
split options from the enlarged pool optimized to gain zero.

`split43_scan.py` also solves body-dependent head polytopes after deleting an
individual medium point. A 300-body pilot and a 3,000-body continuation sampled
all 60 medium families. The latter retained 52,419 exactly checked candidate heads
in 320.28 seconds and produced two two-copy options, both on family 59, row 376,
with the same two tiny deletions. They cannot be selected together and give no
net gain. Their exact options and search receipts are retained.

A useful next constrained search would jointly change the large-family head
configuration or its discriminant orientation before regenerating medium/tiny
splits. Repeating the same fixed large-family polytope with the same templates
has not supplied evidence of a further improvement.

## Reproducibility and storage

All commands must respect the repository CPU reservation rule. This lane used
CPUs 6–7 with one BLAS/OpenMP thread. Generated compatibility inputs, binaries,
large proposal pools and caches are ignored by Git. Final certificates, receipts,
search scripts, short logs and this account of unsuccessful branches are kept.
No shared registry or old verification package was modified by this lane.
The new canonical package contains its own independent verification files.
