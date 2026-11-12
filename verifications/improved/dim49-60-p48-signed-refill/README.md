# Signed P48 classes: exact continuations in dimensions 49–60

The new mother class has **15,622 signed vectors**, of which **10 do not have
an antipode in the class**. Its predecessor had 15,598 signed vectors, arranged
in 7,799 antipodal pairs. Every vector is an integer vector of norm six in the
pinned P48p Gram lattice, and every distinct pair has raw inner product at most
two. Independent bounded-binary64 and scalar-int64 replays check all
**122,015,631** unordered pairs.

This gives **52,431,624** kissing points in 49 dimensions, a gain of 24 over
the preceding certificate. Applying its exchanges to the other direction
classes also produces improvements in 50–60D. These gains are modest; the
structural advance is that antipodality of the mother classes is no longer
required. The published P48p minimum and minimal-shell population remain
mathematical inputs. No unrestricted optimality claim is made.

## Registered bounds

| Dimension | Previous repository bound | New lower bound | Gain |
|---|---:|---:|---:|
| 49 | 52,431,600 | **52,431,624** | +24 |
| 50 | 52,462,800 | **52,462,872** | +72 |
| 51 | 52,509,596 | **52,509,740** | +144 |
| 52 | 52,603,116 | **52,603,402** | +286 |
| 53 | 52,727,584 | **52,728,054** | +470 |
| 54 | 52,975,804 | **52,976,623** | +819 |
| 55 | 53,392,162 | **53,393,546** | +1,384 |
| 56 | 54,260,992 | **54,263,398** | +2,406 |
| 57 | 54,757,822 | **54,760,693** | +2,871 |
| 58 | 56,266,316 | **56,270,076** | +3,760 |
| 59 | 56,948,468 | **56,952,414** | +3,946 |
| 60 | 58,627,711 | **58,631,571** | +3,860 |

## Reproduction

`final.py` exposes `results()` and derives the twelve registered totals from the
pinned parent masks and signed delta sizes. `verify.py` is the complete read-only
entry point: it runs both single-class arithmetic paths, all final-family checks,
and the tail/pole proof, with manifest guards before and after.

From the repository root, on a machine with NumPy and a C++17 compiler:

```sh
OPENBLAS_NUM_THREADS=1 python verifications/improved/dim49-60-p48-signed-refill/verify_single_class.py --cpp
OPENBLAS_NUM_THREADS=1 python verifications/improved/dim49-60-p48-signed-refill/verify_signed_family.py --tails
OPENBLAS_NUM_THREADS=1 python verifications/improved/dim49-60-p48-signed-refill/negative_controls.py
```

The first command independently checks the entire new 49D mother class twice.
The second reconstructs every final class from its predecessor and explicit
signed deltas, freshly checks the complete base geometry, every exceptional
pair, integral isometries, all tail and pole codes, and exact global ownership.
The third rejects a duplicated vector, a wrong norm, and a sign reversal that
preserves norms and distinctness but creates an actual inner product three.
Neither verification command performs an optimization or modifies its inputs.
The JSON receipts are evidence of completed runs; the checkers derive the
counts from certificate coordinates.

The signed deltas are deterministic `.json.gz` files.
`compression_verification.json` proves their decoded bytes are identical to the
JSON inputs used in the completed full pair replays. `certificate.json` pins both
compressed and decoded hashes.

The files in `verifications/improved/dim49-60-p48-pointwise-refill/` supply the
preceding classes, P48 Gram, direction codes and pole codes. Each receipt
records the exact source and delta hashes. No external file, search cache,
network access or expanded 52-million-point mother shell is needed.

## Why a signed class is sufficient

Let M be the complete unit minimal shell of P48p, so |M| = 52,416,000 and
u·v ≤ 1/2 for distinct u,v in M. Choose pairwise disjoint **signed subsets**
C_i of M with u·v ≤ 1/3 for distinct u,v within each C_i. The subsets need
not contain the antipodes of their members.

In a verified antipodal k-dimensional direction code, assign a distinct
antipodal pair {z_i,−z_i} to each C_i. Remove every chosen signed mother
point u from the equator and replace it by

\[
  (\sqrt3\,u/2,z_i/2),\qquad(\sqrt3\,u/2,-z_i/2).
\]

Retain all other mother points at the equator and add a separately verified
unit pole code of cardinality P in the k-dimensional tail. For
S = Σ_i |C_i|, the result has

\[
  |M|-S+2S+P=52,416,000+S+P
\]

points in dimension 48+k. Every pair type has the required angle:

- Within one class, distinct heads give at most (3/4)(1/3)+1/4 = 1/2.
- The two lifts of the same signed head have product 3/4−1/4 = 1/2.
- Across classes, global signed ownership makes the mother heads distinct,
  and distinct antipodal direction pairs have all cross products at most
  1/2. The resulting bound is 3/8+1/8 = 1/2.
- An equator point and a lift have product at most √3/4 < 1/2.
- A lift and a pole have product at most 1/2 by Cauchy–Schwarz.
- Equator/equator and pole/pole pairs inherit their verified angle bounds.

The different tail norms separate equator points, lifts and poles. The
head sets' distinctness and signed ownership prevent duplicate lifted points.
Antipodes may occur in different classes: they are distinct mother points and
satisfy the same proof.

## What is checked exactly

The P48 Gram file is pinned by SHA-256
`01a5d107e4424092e5300e201592aaff62804ee27eeb1d23b6e516249897d3e1`.
The checker confirms all integer lattice coefficients have norm six and
checks each integral transformation T by TGTᵀ = G. The complete old
7,783-line base is checked at absolute product two. This handles every
retained base/base pair, even after individual signs have been removed.
Every other final class vector is checked against all retained signs of
the base and every other exceptional vector.

Global ownership uses complete 48-coordinate integer records, not numerical
hashes. The old classes have disjoint unoriented lines. All removed and added
signed points are explicitly transformed to global coordinates; an added point
cannot be retained elsewhere unless its old occurrence was removed. The
reported total is recomputed from these sets and independently from final
class sizes.

Binary64 products are used only after proving that their inputs, products
and every possible partial sum are integers with absolute value below 2^53.
Consequently these are exact integer calculations with no angular tolerance.
The independent C++ replay uses bounded scalar int64 accumulation and checks
the complete signed class without relying on a conflict graph or search score.

## Discovery

[`signed_p48.py`](../../../research/dim32plus_20261113/middim/signed_p48.py) begins with the preceding 15,598-point symmetric class and
both signs of the retained low-conflict candidate pool. It makes one-for-two
signed exchanges. [`p48_lns.py`](../../../research/dim32plus_20261113/middim/p48_lns.py) then searches complete small blocker
neighborhoods and finds two-for-three exchanges. The final template removes
32 old signed points and adds 56. Its finite candidate pool and termination
are not an exhaustion proof.

[`p48_propagate.py`](../../../research/dim32plus_20261113/middim/p48_propagate.py) applies that explicit 32-for-56 template to each original
class. It tests each new point against all retained points, resolves actual
signed ownership, and retains only positive exchanges. It never multiplies
the single-class gain by the number of classes without checking collisions
and differing original class contents.
