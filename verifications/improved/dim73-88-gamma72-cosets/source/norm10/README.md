# Norm-ten layers over Gamma72

The parent is the even unimodular Gamma72 lattice of minimum squared norm 8.
Its minimal shell has 6,218,175,600 points. These are literature inputs from
[Nebe's lattice catalogue](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/neb72.html)
and [the construction paper](https://arxiv.org/abs/1008.2862).
Every nonempty shell is a spherical 11-design, by Venkov's theorem; see
[Nebe's survey, Theorem 4.1](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/venklat.pdf).

## Geometry

If y has squared norm 10 and u is minimal, minimum 8 applied to y+u and y-u
gives |y.u| <= 5. Thus (4/5)y is compatible with the entire equator. No
Voronoi nearestness condition is needed: the requirement is an inner-product
inequality against the minimal shell.

More generally, for a tail z of squared norm q and a rational number s, use

    p(y,z) = (sqrt(s) y, sqrt((8-10s)/q) z).

Every point has squared norm 8. Equatorial compatibility follows when
s <= 16/25. For distinct head lines, minimum 8 gives |y.y'| <= 6.
At s >= 4/7, head products of absolute value at most three are always safe.
The remaining conflicts between antipodal pairs are determined exactly by

    |s (y.y') + (8-10s) (z.z')/q| > 4.

At s >= 3/5, a pure-tail kissing configuration of squared norm 8 is also
compatible, by Cauchy--Schwarz. We include the 40 D5, 72 E6 or 126 E7 roots.
Their norms, pair products and joint tail-space ranks are checked directly.

## Explicit 73D addition

`shell10_selected.npy` contains 6,235 oriented integral norm-ten vectors y.
For distinct rows, -6 <= y.G.y' <= 3. Against all 6,272 completed half-coset
heads x and all 1,731 old cap lines u, they satisfy |y.G.x| <= 5 and
|y.G.u| <= 3. `verify_shell10.py` checks these statements in int64 arithmetic,
covering 88,773,930 products including diagonal entries.

Add the 12,470 points +/-(4y/5, sqrt(8/5)) to the preceding 73D construction.
Products between upper points lie in [-56/25,88/25], so antipodal completion
is safe. Against the old coset layer the bound is (2/5)*5+2=4. Against an
old cap it is (6 sqrt(3)+4 sqrt(5))/5 < 4; the verifier proves the latter by
integer squarings. Their four absolute latitudes are distinct, so no points
coincide. The old non-lattice witness is retained.

The total is **6,218,204,076**, exceeding the completed local 73D lattice's
6,218,188,144 contacts by **15,932**. The old surgery certificate and completed
lattice certificate are replayed before this addition by `verify73.py`.

## Exact whole-shell existence bounds in 77--79D

These bounds prove the existence of finite configurations. They do not ship
lists of every selected norm-ten head. Only the small tail sets are enumerated.
`norm10_bound_certificate.py` uses exact rational arithmetic and int64 tail
products; it does not rely on the discovery optimizer or floating-point angles.

The theta series is

    E4^9 - 2160 E4^6 Delta + 965520 E4^3 Delta^2 - 27302400 Delta^3.

The checker derives these coefficients by cancelling q, q^2 and q^3, then
computes N8=6,218,175,600 and N10=15,281,788,354,560. For the norm-r shell
and a squared-norm-t vector x, the 11-design identities give

    sum (v.x)^(2j) = Nr (2j-1)!! (rt)^j / (72*74*...*(72+2j-2)),  0<=j<=5.

Applying this to the minimal shell against a norm-ten vector determines
the absolute-product populations at 0,...,5. The population at 4 is 2,112,250.
The map x -> y-x bijects norm-ten vectors with y.x=6 and minimal vectors
with y.(y-x)=4. Thus 2,112,250 is also the absolute-product-six population
of the norm-ten shell. Its six remaining populations follow from the same
moment identities, after subtracting the two vectors +/-y:

| Absolute product | Population |
|---:|---:|
| 0 | 5,118,989,848,896 |
| 1 | 7,238,164,864,000 |
| 2 | 2,503,639,594,500 |
| 3 | 395,496,128,000 |
| 4 | 24,994,840,000 |
| 5 | 500,966,912 |
| 6 | 2,112,250 |
| 10 | 2 |

Choose one orientation per head line. Graph vertices are pairs (head line,
tail), representing the antipodal pair +/-p(y,z). Two vertices are adjacent
exactly when their pairs conflict. All supplied tails are antipodal. Therefore
the number of conflicts against a head line with absolute product j is
independent of its chosen orientation. If b_j(z) counts tails w satisfying

    (8-10s) z.w > (4-js) q,

then the closed degree is

    b_10(z) + a4*b_4(z)/2 + a5*b_5(z)/2 + a6*b_6(z)/2.

Here b_10 is the closed contribution from the same head line. For the three
reported constructions it equals the entire tail population. There are a_j/2
other head lines of absolute product j. All divisions by two are exact.

A uniformly random ordering selects a vertex when it is first in its closed
neighborhood. Selected vertices form an independent set, and each vertex is
selected with probability the reciprocal of its closed degree. Taking
expectations and rounding up proves the Caro--Wei bound used by the checker.
The number of points added is twice this integer bound, plus the pure tails.

| Dimension | s | q | Tails | Earlier guaranteed pairs | Pure tails | Current weighted total |
|---:|:---:|---:|---:|---:|---:|---:|
| 77 | 16/25 | 40 | 560 | 144,688 | 40 | **6,218,490,542** |
| 78 | 16/25 | 40 | 2,160 | 255,104 | 72 | **6,218,688,770** |
| 79 | 12/19 | 32 | 4,158 | 368,155 | 126 | **6,218,922,650** |

The tails are sections of 2E8: in rank five the last three coordinates vanish;
in rank six the first three are equal; in rank seven the first two are equal.
Completeness of these tail shells is unnecessary. The checker verifies the
supplied finite sets, their exact ranks, antipodality, profiles and the rational
bound. `norm10_tails.npz` and `norm10_bound_certificate.json` hold the data and
the replay report. The package wrapper `verify_norm10.py` checks both the
73D addition and these whole-shell proofs in a temporary copy.

## Mixed tails and unequal priorities, 28 October 2026

`weighted_norm10_certificate.py` strengthens the same existence proof. It
retains the full Gamma72 equator and the pure-tail roots, mixes compatible
candidate families, and replaces uniform random priorities by fixed integer
rates. Run `python ../../verify_weighted_norm10.py` from this directory for
the isolated replay, or `python weighted_norm10_certificate.py --check` to
recompute the saved report without writing it. All 22,848,388 ordered tail
products are checked; diagonal entries are included in that count.
`verify_weighted_independent.py` reconstructs the tail sets from scratch,
checks all design-moment identities, and recomputes the neighborhoods
with integer square-root cutoffs instead of squared masks. It also
enumerates all 64 graphs on four vertices and all 24 weighted clock
orders, verifying the local-minimum expectation and independence.
Both isolated wrappers run this independent check as well.

| Dimension | Tail norms q | Tails | Guaranteed pairs | Pure tails | Current bound | Gain over preceding bound |
|---:|:---:|---:|---:|---:|---:|---:|
| 77 | 32, 40 | 760 | 157,451 | 40 | **6,218,490,542** | 25,526 |
| 78 | 8, 40 | 2,232 | 256,549 | 72 | **6,218,688,770** | 2,890 |
| 79 | 32 | 4,158 | 373,462 | 126 | **6,218,922,650** | 10,614 |

The rank-five tails add all 200 vectors `(2v,0,0,0)` with integral `v` in
five coordinates, `v.v=8` and even coordinate sum. Their 40 two-coordinate
vectors receive rate 1050 and the remaining 160 receive rate 818. Among
the preceding norm-40 tails, the 80 vectors having a coordinate of absolute
value six receive rate 882 and the other 480 receive rate 1000.

The rank-six tails add the 72 pure-tail roots as *tail directions* at rate
1015, while the preceding 2160 norm-40 tails receive rate 1000. These
candidate antipodal pairs have nonzero head components and differ from
the pure-tail root points, which are also retained.

In rank seven, the 126 norm-32 tails equal to twice a pure-tail root receive
rate 1000; the other 4032 tails receive rate 311. The head scales remain
16/25, 16/25 and 12/19 in ranks five, six and seven. Rates are invariant
under negating a tail. No optimizer, floating-point threshold, orbit
transitivity assumption, or new external coordinate file is trusted by
the replay.

### Exact mixed-shell conflict test

Normalize each tail separately: `t_z=z/sqrt(q_z)`. For different head lines
of absolute product `j`, a positive head product conflicts exactly when

    (8-10s) (z.w)/sqrt(q_z q_w) > 4-js.

The thresholds are positive for `j=4,5,6` in these constructions. Writing
`s=p/d`, this is equivalent to both integer inequalities

    z.w > 0,
    (8d-10p)^2 (z.w)^2 > (4d-jp)^2 q_z q_w.

Equality is compatible. Products at absolute head value at most three
are harmless since `3s+(8-10s)<=4`. Every tail choice above a *single*
head line conflicts with every other such choice: even the smallest
product is `10s-(8-10s)>4`. Its closed-neighborhood contribution therefore
contains every tail, including itself.

Let `r(z)` be the supplied positive integer rate and put

    B_j(z) = sum_w r(w) [the strict mixed-shell conflict test holds],
    D(z) = sum_w r(w) + (a4 B_4(z)+a5 B_5(z)+a6 B_6(z))/2.

Negating all tails preserves rates, so `B_j` is also the conflict rate
when the chosen orientations of two head lines have negative product.
There are `a_j/2` other head lines of each absolute product. Thus `D(z)`
is the exact sum of rates in the closed neighborhood of `(head line,z)`.
The checker computes it for **every** tail, without extrapolating from
representatives. Tail and pure-root spans are checked together by an
exact rational rank computation. The earlier equatorial and pure-tail
Cauchy--Schwarz bounds apply to all normalized tails.

### Weighted-priority proof

Give each graph vertex an independent exponential clock with its assigned
positive rate. Keep a vertex exactly when its clock is smaller than all
its neighbors' clocks. Two adjacent vertices cannot both survive. For
independent exponential clocks, the probability that a given clock rings
first in a finite set is its rate divided by the sum of rates in that
set. Consequently the expected size of the resulting independent set is

    alpha = (N10/2) sum_z r(z)/D(z).

At least one ordering has at least `ceil(alpha)` vertices. Multiplying
by two and adjoining the equator and pure tails proves the displayed
bounds. The report stores each rational expectation, every `(rate,D)`
profile and its multiplicity; its ceiling is computed with integers.
This certifies existence without constructing the enormous head set.
It does not assert optimal rates, optimal tails, or a construction ceiling.

### Relation to PackingStar

The [published PackingStar method](https://arxiv.org/abs/2511.13391) searches
Gram matrices over restricted cosine sets, using a filler to add candidates
and a corrector to remove poor choices while respecting positive
semidefiniteness and dimension. Its earlier Leech construction also found
better zero-sum triangle partitions. The
[26 October repository update](https://github.com/CDM1619/PackingStar/blob/67be454b83af4287e1a2b61e5bd8c3596a8061d6/README.md)
describes joint optimization of removals, midpoint insertions and auxiliary
layers, with cosines in multiples of 1/8 or 1/16. The public snapshot
contains configurations and verifiers, but no training/search implementation
for this newer filling scheme.

The useful transferable principle is to vary candidate families and their
selection jointly. Here that leads to mixed normalized tail shells and
nonuniform conflict selection, which the old single-shell, uniform-order
proof omitted. A dense Gram matrix for over six billion contacts is
unnecessary: the design identities describe the head shell symbolically,
and only the small tail Gram matrices are enumerated. This is an independent
exact adaptation of that search principle; it does not reproduce their
reinforcement-learning system. Improvements are measured against this
repository's preceding certificates. No new global priority survey is claimed.

## Further progress in 73D and 82--87D, 28 October 2026

`refilled_shell10.npz` replaces the 6,235-head explicit 73D layer with
**7,663** oriented norm-ten heads. Its 15,326 points give **6,218,206,932**,
an increase of **2,856**. The equator, half-coset contacts, cap surgery and
non-lattice witness remain. Integral automorphism images, two-for-one
replacements and neutral exchanges supplied the candidates; the final
coordinates alone are trusted. `verify_refilled_shell10.py` checks all
**120,048,558** head/head and mixed products in int64 and independently
in binary32 after proving every term and partial sum is an exactly
represented integer below 2^24. `refilled_shell10_verification.json` and
`refill_search.json` record the checks and finite search scope.

The signed-tail certificates improve five more dimensions:

| Dimension | Rank k | Support weights | Guaranteed pairs | Pure tails | Current bound | Gain over preceding repository bound |
|---:|---:|:---:|---:|---:|---:|---:|---:|
| 82 | 10 | 1, 4 | 997,443 | 180 | **6,220,170,666** | 169,130 |
| 83 | 11 | 1, 4 | 1,416,419 | 220 | **6,221,008,658** | 474,298 |
| 84 | 12 | 1, 4 | 1,924,770 | 264 | **6,222,025,404** | 468,208 |
| 85 | 13 | 1, 4, 9 | 2,523,388 | 312 | **6,223,222,688** | 1,217,848 |
| 87 | 15 | 1, 4, 9 | 4,009,295 | 420 | **6,226,194,610** | 356,580 |

Let T_w contain all vectors in R^k with exactly w coordinates in {+1,-1}
and every other coordinate zero. It has 2^w binomial(k,w) vectors of
squared norm w. Take the indicated union of T_1, T_4 and T_9, and use

    +/-(sqrt(8/13) y, sqrt(24/(13w)) z), y.y=10, z in T_w.

These points have squared norm eight. Equatorial products are at most
5 sqrt(8/13)<4. Head products of absolute value at most three give at
most 48/13<4. All choices above one head line conflict, since even the
smallest positive-body product is 56/13>4. Adjoin all 2k(k-1) pure roots
`2(+/-e_i +/-e_j)`, i<j. Their norm is eight, distinct products are at
most four, and products against the new layer are at most sqrt(192/13)<4.
T_1 spans R^k, and every other tail and pure root lies in that same space.

Against a fixed positive weight-a tail, a weight-b tail overlaps it in r
positions with l positive signs, giving product d=2l-r and population

    binomial(a,r) binomial(k-a,b-r) binomial(r,l) 2^(b-r).

Sum over valid r,l to obtain H_ab(d). Signed coordinate permutations
act transitively on T_a, preserve each T_b and its assigned rate, and
preserve products. Thus this is the exact profile for every tail in T_a.
For absolute head products j=4,5,6, define the strict conflict populations

    B_j(a,b) = sum_d H_ab(d) [d>0 and 24^2 d^2 > (52-8j)^2 ab],
    M_ab = |T_b| + sum_(j=4,5,6) (a_j/2) B_j(a,b),
    D_a = sum_b M_ab rate_b,
    alpha = (N10/2) sum_a |T_a| rate_a/D_a.

Boundary products are compatible. Negating tails preserves rates, so
negative head products have the same populations. M_ab includes the full
same-head closed neighborhood. The preceding exponential-clock argument
therefore guarantees ceil(alpha) antipodal pairs. Use the integer rates:

| Rank k | T_1 rate | T_4 rate | T_9 rate |
|---:|---:|---:|---:|
| 10 | 19,871,466 | 10,000,000 | unused |
| 11 | 21,748,661 | 10,000,000 | unused |
| 12 | 23,745,110 | 10,000,000 | unused |
| 13 | 26,122,963 | 10,000,000 | 3,247 |
| 15 | 31,936,463 | 10,000,000 | 4,877 |

`signed_norm10_certificate.py` recomputes every combinatorial population,
neighborhood matrix, weighted degree and rational expectation; the report
is `signed_norm10_certificate.json`. `verify_signed_profiles.py` instead
computes H_ab from coefficients of

    (1+t(X+X^-1))^a (1+2t)^(k-a), at t^b X^d,

using integer dynamic programming and integer square-root cutoffs. Direct
enumeration of every signed orbit in ranks one through five checks this
polynomial model. The independent checker also verifies every pure-root
product and exact rank. No large tail list, selected head list, optimizer
output or floating-point threshold is a proof input. The Gamma72 minimum,
theta series and all-shell 11-design facts remain the literature inputs.
Run `python ../../verify_norm10.py` for the complete isolated replay.
