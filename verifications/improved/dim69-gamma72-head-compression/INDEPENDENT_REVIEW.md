# Adversarial independent review — 13 November 2026

**Verdict: the claimed lower bound τ(69) ≥ 743,522,304 passes this review.**
The gain over the preceding registered 689,452,922 is 54,069,382 points.
No mathematical flaw was found in the finite proof. The published Gamma72
catalogue identity, minimum 8, minimal-shell population 6,218,175,600, and
spherical 11-design property remain explicit inputs.

## Independent arithmetic reconstruction

The third review traversed the normal-coordinate generators on individual
profiles, instead of using the producer's canonical orbit representatives.
It recovered all 1,681 permitted profiles, 24 known normal roots, and 16
averaging orbits. Every generator preserves both the root set and the profile
domain. Because the roots span the normal space, their finite permutation
action is faithful, which justifies finite orientation averaging.

Sphere moments were independently recomputed using factorial expressions,
and all ten primitive integer moment rows and right-hand sides were reproduced.
The root-basis and orthonormal profile descriptions, heights, and squared
rescalings agree exactly. Rational weak duality gives **743,522,304 exactly**.
Every one of the 55 geometric pair types passes. The smallest strictly
positive squared margin is **1/180**.

## Points examined adversarially

**Averaging without an ambient symmetry.** The argument averages actual normal
tuples related by automorphisms of the embedded 2D4 sublattice. Each tuple is
formed from integer combinations of existing lattice vectors. It does not
assert that these transformations preserve the entire Gamma72 lattice or
its shell. The largest actual ten-slice population is at least their average.

**Completeness of the profile domain.** Integral products with the roots
2(e_i±e_j) imply that h=4p has integer coordinates of common parity.
Every minimal vector other than one of the 24 normal roots has absolute
product at most 4 with each normal root. These are exactly the stated
top-two-coordinate inequalities. The 24 exceptional roots are separately
subtracted from every moment equation. A source point cannot have another
unaccounted normal-root profile: that profile already has squared norm 8,
so its orthogonal component would have to vanish.

**Projection and image collisions.** On a fixed profile, the orthogonal
component determines the original vector uniquely, and the rescaling is
positive. Across different profiles, the exact pair inequalities apply to
distinct source points. A coincident image would have product 8 with itself,
contradicting the proved upper bound 4. The population therefore survives
the construction without hidden duplicate removal.

**Signs before squaring.** The selected profiles make both sides of the
unsquared geometric comparison nonnegative. The rational squared tests are
therefore sufficient. The review checked these signs independently.

## Semantic corruption controls

Run `python rejection_controls.py` in this directory. It uses only
package-local files and the Python standard library, creates temporary copies,
and leaves the mathematical inputs unchanged. Both unmodified certificates
pass and **26 corrupted or unsafe-invocation cases are rejected**.

The tests cover false counts and populations, duplicate slices, invalid
profile domains, omitted cells, false moments and objectives, altered dual
coefficients, invalid scales, and geometric incompatibility. One test changes
the first height to −1/2 and its squared scale to 1. All within-slice image
norms remain valid, but its cross-slice squared margin with slice 5 becomes
**−31/16**; both checkers reject it.

The altered-normal tests deliberately refresh the embedding hash in an
isolated copy of the checker, so rejection reaches the integer Gram calculation
rather than merely a checksum mismatch. Both checkers also reject execution
with Python's `-O` flag. The review initially found that the primary checker's
`assert __debug__` guard would itself disappear under `-O`; it was replaced
with an explicit `if`/`raise`, and the final controls confirm the fix.

`independent_review.json` records the independent arithmetic findings and
the hashes of the final checked inputs. `rejection_controls.json` lists each
control and its observed result. These checks establish the finite application
of the cited lattice premises; they do not reprove the published lattice
theorems or enumerate the 743,522,304-point image.
