# Structured 21D attack, 4 November 2026

The verified lower bound remains **31560**. This continuation replaces the
previous poorly calibrated point-level search with a smaller hole-set model,
explicit recovery benchmarks, nonlinear starting configurations, and exact
local exclusions. No global impossibility result is claimed.

The construction and comparison bound are due to Kevin Li,
[Improved Kissing Numbers in Nineteen and Twenty-One Dimensions](https://hexagonmath.org/2610.00248v1),
Section 3, 14 October 2026, CC BY 4.0. Inputs are the independent reconstruction
in `../compute19_21_20261104`. Our nonlinear examples are not claimed to be
new up to geometric equivalence.

## The mathematical model

Keep the 27720-point odd-sign base. The auxiliary family has 4096 punctured
Golay sign words and 21 line profiles, with magnitude 2 on a line and 3
elsewhere. Its compatibility rule is exact: if selected words differ by a
line P, their two labels and P must be distinct. The word graph splits into
four isomorphic 1024-word components. A component with 961 selected words
would give 31564 points after replication; the known construction has 960.

Let H be the omitted words, called holes, in one component. Every selected
word u labeled P requires u+P in H. Therefore H is a dominating set in the
21-regular word graph. This is a **necessary** condition for every labeling.

If H is also independent, domination is **sufficient**. Label each u outside
H by any direction P leading to a hole. If adjacent selected u and v both
had label P, their target holes u+P and v+P would be adjacent, contradicting
independence. Neither selected endpoint can point along their joining edge
because the target would not be a hole. Thus all three labels are distinct.

This sufficient subfamily needs 1024 hole indicators instead of 21504
word-label candidates. It is broader than linear kernels, but it is still
a restriction of the full line-label model. `solve_holes.py --relax` drops
independence and searches a necessary relaxation. `solve_sat.py --mode labels`
encodes the full rule, without assuming independence or translation-periodic
labels. All models leave the geometric base and profile family fixed.

## Recovery and verification

`recover.py` samples four binary linear functionals until the 21 line images
cover all 15 nonzero syndromes and avoid zero. It never uses the published
map or labels. Their common kernel has 64 holes, and the theorem above gives
960 selected words. This is a structure-aware recovery method, not a blind
general-purpose optimization benchmark.

- All **32/32** independent-seed recovery tests passed; the slowest took
  **0.147 seconds**, excluding the subsequent pair check.
- Every recovered component was checked by direct integer inner products.
- A generic CP solver with no hints still failed to recover 64 holes in
  120 seconds. The reduction alone does not cure generic search weakness.
- An independent SAT encoding recovered the 64-hole configuration in
  19.44 seconds with phase hints. This is explicitly a warm-start control.
- `quotients.py` exhaustively checked all **651** four-dimensional translation
  subgroups inside the first recovered kernel. In each 64-word quotient it
  enumerated independent dominating sets of at most four holes containing
  zero. Translation makes that last normalization harmless. The minimum was
  four, and the enumeration found **eight distinct nonlinear 64-hole sets**
  after lifting. The 2758 quotient solutions count repeats across periods.
- `verify_geometry.py` replicated one nonlinear component, exported all
  **31560 integer rays**, and checked all **498001020 unordered pairs**.
  There were no violations. Deliberately invalid labels and a duplicated
  ray were rejected. `data/nonlinear21_verified.npz` is the complete artifact.

## A proved barrier to small changes

`certify_trade_barrier.py` checks the following exact statistics for all 40
linear and nonlinear seeds. Each old hole has 12 private neighbors, meaning
vertices covered by that hole and by no other old hole. A possible new hole
has at most three old-hole neighbors, at most 17 private neighbors in total,
and at most two private neighbors belonging to any one old hole.

Suppose a smaller dominating set removes a and b old holes from the two
bipartition classes, and adds c and d new holes. Necessarily

```
c+d <= a+b-1,
12a <= min(2ac,17c)+d,
12b <= min(2bd,17d)+c,
a <= 3d,    b <= 3c.
```

For example, the 12a private vertices left by removed holes in the first
class must be covered by new holes. The c new holes in that same class cover
at most two per removed hole and at most 17 each in total. The d new holes
in the other class can cover at most d of these vertices by being those
vertices themselves. The final two inequalities count domination of the
removed holes themselves, which had no old-hole neighbors.

Exhausting these small integer inequalities gives **no possibilities with
a+b <= 12**. Hence every improvement around any of these 40 configurations
must remove at least **13 existing holes**, even if the new holes are
adjacent. This is a statement about every nearby full line labeling, not
just independent hole sets. It allows all labels to change.

This certificate uses integer counts and elementary inequalities, with no
SAT/LP solver or floating-point tolerance. It explains a real limitation
of small local repairs; it is not a global ceiling for the family.

## Rational certificates for retaining half a kernel

`certify_halves.py` examines the 63 five-dimensional subgroups T of the first
recovered six-dimensional kernel. If all 32 words in T remain holes, average
any prospective hole indicator over T. This yields 32 fractional coset
indicators; it does **not** impose periodicity on the prospective solution.

The quotient has an integer closed-neighborhood matrix A. Its zero-coset
indicator is one. The saved certificate gives nonnegative rational weights
y on rows with A[i,0]=0, with sum(y)=1 and every other column load at most
one. Weighting the domination inequalities forces at least one additional
unit of coset hole mass, hence at least 64 holes in the full component.

Certificates were found for **62 of the 63 subgroups**. Translation also
certifies their other cosets: 124 halves inside the original kernel, or
1984 fixed cosets in the entire component. The exact fraction checker
replays every certificate without an LP solver. `--half-cuts` uses the 124
halves inside the original kernel as necessary exclusions when seeking a gain.

The remaining subgroup, mask 60 in the saved basis, was investigated
separately. Its domination LP bound is approximately 62.439, so that
relaxation permits 63 holes. CP-SAT ruled out 63 holes when this subgroup is
retained **and the new hole set is independent**, in 0.312 seconds. This is
a solver-reported restricted exclusion, not a separately checked rational
or SAT proof. Adjacent-hole solutions require the full label test.

## Searches and interpretation

The run-level records are `data/*_result.json`, with full CP logs alongside.
They include unrestricted independent-domination and domination searches,
full-label SAT searches, radius-8/13/16 changes, nonlinear initial conditions,
the exceptional retained half, and global optimization from a nonlinear seed.
Exact private-neighbor inequalities and rational half-kernel exclusions are
fed back into the later CP models. The global optimization uses a verified
64-hole incumbent and allows every hole indicator to change.

The final 300-second global run retained 64 holes and returned **FEASIBLE**,
with a reported lower bound of 49 holes in the independent-domination model.
The unrestricted full-label SAT run and the exceptional-half full-label run
both timed out after 300 seconds without a solution or an infeasibility proof.
The larger radius-13/16 neighborhoods also remain unresolved. No 961-point
component was found, and these timeouts provide no global ceiling.

`UNKNOWN` and `FEASIBLE` with a remaining gap are timeouts, not exclusions.
The final aggregate is `data/summary.json`; it checks that every run log has
a corresponding result before reporting completion. Valid saved solutions
are checked again against the exact line rule. Formal negative claims here
are limited to the integer trade certificate, rational conditional bounds,
and the stated finite quotient enumeration. CP infeasibility results are
identified separately.

The useful change in practice is that recovery now follows the algebraic
construction, and the search knows some barriers it must cross. Neither the
number of local moves nor these restricted certificates establishes that
31560 is optimal. The target of 961 selected words remains the discovery test.

## Replay

Use Python with NumPy, SciPy, OR-Tools 9.15 and python-sat 1.9.dev20 installed.
The original campaign used `/usr/bin/python3` and installed the latter two in
the enclosing workspace's `.compute-deps`.
Every compute command reserves CPU 127, following the repository AGENTS.md.
This work is discrete and exact-arithmetic focused; no GPU was needed.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
taskset -c 0-126 python recover.py
taskset -c 0-126 python quotients.py
taskset -c 0-126 python certify_trade_barrier.py
taskset -c 0-126 python certify_halves.py --verify-only
taskset -c 0-126 env OPENBLAS_NUM_THREADS=8 python verify_geometry.py
taskset -c 0-126 python solve_holes.py --name global_replay \
  --cap 64 --seconds 300 --workers 24 --hint --optimize --nonlinear \
  --trade-cuts --half-cuts
taskset -c 0-126 python solve_sat.py --name labels_replay \
  --mode labels --hint --seconds 300
taskset -c 0-126 python audit.py
```

Replaying parallel optimization is not guaranteed to reproduce identical
search trajectories. Certificate checks and geometric admission do not
depend on those trajectories. Main published-bound registries were not changed.
