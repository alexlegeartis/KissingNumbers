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

| Dimension | s | q | Tails | Guaranteed pairs | Pure tails | Total |
|---:|:---:|---:|---:|---:|---:|---:|
| 77 | 16/25 | 40 | 560 | 144,688 | 40 | **6,218,465,016** |
| 78 | 16/25 | 40 | 2,160 | 255,104 | 72 | **6,218,685,880** |
| 79 | 12/19 | 32 | 4,158 | 368,155 | 126 | **6,218,912,036** |

The tails are sections of 2E8: in rank five the last three coordinates vanish;
in rank six the first three are equal; in rank seven the first two are equal.
Completeness of these tail shells is unnecessary. The checker verifies the
supplied finite sets, their exact ranks, antipodality, profiles and the rational
bound. `norm10_tails.npz` and `norm10_bound_certificate.json` hold the data and
the replay report. The package wrapper `verify_norm10.py` checks both the
73D addition and these whole-shell proofs in a temporary copy.
