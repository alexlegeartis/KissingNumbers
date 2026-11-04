# Certified upper bound for the 21D line-profile family

The full 1024-word component satisfies

**960 <= M <= 971.**

The upper bound is now an exact, independently checked certificate. It
improves the earlier solver-reported upper bound of 975. This work does
**not** prove M=960 and does not produce a configuration with more than 960
points in a component.

For the fixed 27720-point base and all four components, the corresponding
family interval is **31560 to 31604**. The upper endpoint is a restriction
of this construction family, not an upper bound on the 21D kissing number.

The geometry and known construction are Kevin Li's
[line-weighted Golay construction](https://hexagonmath.org/2610.00248v1),
Section 3, CC BY 4.0. The input geometry was independently reconstructed in
`../compute19_21_20261104`. The general use of symmetry and moment matrices
is related to [Schrijver's semidefinite bounds for codes](https://homepages.cwi.nl/~lex/files/codes.pdf).
The graph-specific models and certificates here were computed in this campaign.

## Exact result and verifier

`data/upper971_certificate.json` contains an exclusion of an independent
set of size 972 in the complete word-label conflict graph. It uses 42
inequalities and six equalities with rational multipliers. Their weighted
sum cancels every variable and gives the contradiction

```
0 >= 10710060036218601246681000926260048316204797557775338473473059284269
     / 18601062554535530604433104378115630932512658256839047126543708641243680
  > 0.
```

The verifier `verify_upper.py` requires only Python, NumPy and the saved
inputs. It uses no LP, SDP or SAT solver, numerical eigenvalue test, or
optimization tolerance. It independently checks:

- The geometric conflict rule on all 451584 word-difference/label cases.
- Each supplied graph automorphism, its linearity under XOR, and the
  induced permutation of the 21 labels.
- Connectivity of every claimed pair-moment orbit, ensuring that inequivalent
  pairs were not accidentally merged.
- Every clique and positive-semidefinite quadratic-form inequality used in
  the proof, and exact cancellation of their rational weighted sum.
- Rejection of a corrupted certificate.

The actual known 960-point component also passes all 559 intermediate
inequalities, using exact averaged ordered-pair counts. See
`data/baseline_control.json`.

## Why the certificate covers arbitrary labelings

There are 21504 candidates `(word,label)` in each component. A word has 21
profiles. Different profiles at the same word conflict. If two different
words differ by a line P, their candidates conflict precisely when either
label equals P or the two labels coincide. Other word differences are
compatible. This gives the degree-81 conflict graph.

Suppose an independent set I has exactly 972 vertices. Average its outer
product indicator matrix over word translations and verified graph
automorphisms, and divide by the common selection probability. This matrix
is positive semidefinite, entrywise nonnegative, and has diagonal one.
Every edge entry is zero and every row sums to 972.

This averaging imposes no symmetry on I itself. The group action is
transitive on vertices. After translation, entries depend on a word
difference d and two labels a,b. The verified automorphisms and matrix
symmetry give 48 entry classes, denoted f_j.

For every binary character chi of the 1024-word translation group, the
21-by-21 Fourier block

```
B_chi[a,b] = sum_d (-1)^(chi dot d) f[d,a,b]
```

is positive semidefinite. Consequently every **integer** vector v gives
the valid linear inequality `v^T B_chi v >= 0`. The proof saves each v and
chi explicitly; its coefficients are rebuilt using integer arithmetic.
Numerical eigensolvers were used to discover useful vectors, but their
accuracy is irrelevant to the validity of the resulting inequalities.

Another valid inequality comes from a 22-vertex clique: all 21 candidates
at a word d, plus the candidate `(d+P,P)`. Multiplying its independent-set
inequality by the indicator of a fixed candidate and averaging gives a
linear inequality in the pair moments with right-hand side one. These
clique inequalities, nonnegativity, the Fourier inequalities, edge zeros,
diagonal normalization and the fixed row sum suffice for the certificate.

Any independent set larger than 972 would contain one of size 972. Thus
excluding that size proves M <= 971. The proof includes every label and
does not require independent omissions, periodicity, proximity to a known
configuration, or a fixed subset of omitted words.

## An exact limitation of the tested pair-moment relaxation

`data/relaxation971_rational.json` contains an exact feasible **moment
solution at cardinality 971**. This is not an independent set or a point
configuration. Its role is to demonstrate that the tested relaxation
really does permit 971.

`verify_relaxation.py` checks it using exact fractions and rational LDL
factorization with positive pivots. The verifier covers all 1024 Fourier
characters through eight verified character orbits. It also checks all
21504 conditional 22-clique inequalities, 430080 conditional four-clique
inequalities, and products of the maximal-clique slack indicators.

Therefore this pair-moment SDP, even with all these tested clique
strengthenings, **cannot certify an integer upper bound below 971**.
More iterations or tighter floating-point tolerances cannot make it prove
960. Higher-order information or a discrete exclusion is needed. This
statement concerns the specified relaxation, not all possible SDP methods.

## Other global models tested

The domination route was investigated first. A selection of 961 words
requires 63 omitted words, and those omissions must dominate the word graph.
An independent dominating omission set would suffice for compatible labels,
but independence was not assumed in the main upper-bound models.

The word graph has eight verified vertex orbits under its zero stabilizer,
168 ordered-pair orbits, and 42 unordered-triple moment orbits. We built
degree-three covering inequalities by multiplying domination constraints
by pairs of Boolean literals. Fixed-cardinality equations, Fourier
positivity and conditional moment positivity were added.

`moments.py` implements these global necessary conditions. `algebra.npz`
contains integer multiplication matrices for the invariant matrix algebra;
its Gram representation reduces the two conditional PSD tests from order
1024 to order 168. Numerical results permit 63 holes, including when
independence is additionally imposed. A rational LP contradiction at 49
holes was also extracted (`farkas49.json`), but the standalone full-label
certificate above is stronger and is the result relied on here.

`coupled.py` links the label-pair moments to the hole-pair marginals of this
degree-three model. The tested numerical relaxation remains feasible at
971 selected words. Those coupled numerical outputs are diagnostic only;
they are not claimed as exact feasible certificates.

No new large SAT campaign was counted as progress toward this bound. The
improvement comes from stronger global inequalities and an exact proof,
not another timeout or a count of search attempts.

## Files and replay

The authoritative results are `data/upper_verification.json`,
`data/relaxation_verification.json`, and `data/summary.json`.
`data/manifest.json` hashes source, input and result artifacts.

Certificate replay needs no optimization package:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
taskset -c 0-126 python verify_upper.py
taskset -c 0-126 python verify_relaxation.py
```

Discovery uses NumPy, SciPy, SymPy, CVXPY 1.7.5, Clarabel 0.11.1 and pynauty
2.8.8.1. Install these in the Python environment used below. The original
campaign installed extra packages in the enclosing workspace's `.compute-sdp`
directory. The C++ orbit tools are regenerable from source.

```sh
taskset -c 0-126 g++ -O3 -std=c++17 orbits.cpp -o orbits.exe
taskset -c 0-126 g++ -O3 -std=c++17 label_orbits.cpp -o label_orbits.exe
taskset -c 0-126 python prepare.py
taskset -c 0-126 python labels_sdp.py --cliques --four
taskset -c 0-126 python label_outer.py --exclude 972
taskset -c 0-126 python make_farkas.py
taskset -c 0-126 python export_certificate.py
taskset -c 0-126 python verify_upper.py
taskset -c 0-126 python audit.py
```

The upper-bound outer approximation generates exact integer PSD cuts, then
uses floating-point LP only to discover a useful combination. The saved
rational combination is exported by `export_certificate.py` and verified
independently. All computation reserved logical CPU 127. This campaign
used CPU algebra and optimization; it did not require a GPU. No published
kissing-number registry was modified.
