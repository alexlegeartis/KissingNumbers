# Pointwise P48 refill in dimensions 49–60

**Retained construction.** The twelve bounds below are the unchanged preceding
certificates. The [13 November signed refill](../dim49-60-p48-signed-refill/)
now supplies the current registry values in dimensions 49–60, including
**58,631,571** in 60D. Historical gains below retain their original comparisons.

Twelve explicit constructions improve the repository's preceding lower bounds.
The 60D construction has **58,627,711 points**, an increase of **84,186**.
The finite certificates are antipodal mother-line classes, enlarged by pointwise
insertions and ownership reassignment. They need not remain subsets of the
original 7,783-line images. These are lower bounds, not optimality or worldwide-priority claims.

| Dimension | Previous bound at publication | Retained bound | Historical gain | Selected mother lines | Current registry |
|---|---:|---:|---:|---:|---:|
| 49 | 52,431,568 | **52,431,600** | +32 | 7,799 | 52,431,624 |
| 50 | 52,462,704 | **52,462,800** | +96 | 23,397 | 52,462,872 |
| 51 | 52,509,404 | **52,509,596** | +192 | 46,792 | 52,509,740 |
| 52 | 52,602,708 | **52,603,116** | +408 | 93,546 | 52,603,402 |
| 53 | 52,726,898 | **52,727,584** | +686 | 155,772 | 52,728,054 |
| 54 | 52,974,448 | **52,975,804** | +1,356 | 279,866 | 52,976,623 |
| 55 | 53,389,510 | **53,392,162** | +2,652 | 488,018 | 53,393,546 |
| 56 | 54,254,038 | **54,260,992** | +6,954 | 922,376 | 54,263,398 |
| 57 | 54,747,366 | **54,757,822** | +10,456 | 1,170,758 | 54,760,693 |
| 58 | 56,237,990 | **56,266,316** | +28,326 | 1,924,903 | 56,270,076 |
| 59 | 56,907,964 | **56,948,468** | +40,504 | 2,265,932 | 56,952,414 |
| 60 | 58,543,525 | **58,627,711** | +84,186 | 3,105,435 | 58,631,571 |

## Portable verification

From this directory, run:

```sh
python verify.py
python test_certificate.py
```

Python 3.9+, NumPy and threadpoolctl are required. The full positive replay is
read-only and needs no network, optimizer, outside coordinate files, search
cache or GPU. It checks the manifest before and after the proofs, derives counts
from the final masks, and runs two distinct exact head verifiers. The control
suite uses temporary files and rejects twelve malformed-input or disabled-check
cases. `final.py` is only a count driver, not a geometry proof.

- `verify_heads.py` uses bounded int64 products in each image's original frame.
- `independent_audit.py` first materializes each final class in global coordinates,
  then uses binary64 products whose inputs, products and every partial sum are
  proved to be integers of absolute value below 2^53. These are exact integer
  computations, not tolerance-based numerical tests.
- Both prove all 30,283,653 original base pairs, exact norms and isometries, every
  added/retained and added/added class pair, and complete-coordinate global
  ownership. No probabilistic hashes are used for disjointness.
- `tail_and_lift.py` checks all 751,195 direction-code pairs and all 353,220 pairs
  in the independent 841-point pole code, then checks every lift pair category.
- `audits/` contains the independently generated receipts and rejection controls.
  Receipts are evidence of a run, not substitutes for the checkers.

## Construction

Work in the published P48p Gram lattice. Its minimum squared norm six and
52,416,000 minimal vectors are inherited published facts. A row x of integer
lattice coordinates with x G x^t = 6 specifies a unit mother vector after
normalization in this Gram metric. No fresh enumeration of this shell is claimed.

For class i, take the base rows selected by `retained_bits[i]`, append the rows
in `extra_points[extra_offsets[i]:extra_offsets[i+1]]`, and right-multiply by
`transforms[i]`. Each stored row denotes both signs of its mother line. The
verifiers prove distinct lines within a class have absolute raw product at most
two, and that no mother line occurs in two classes.

In dimension 48+k, assign the classes to distinct antipodal pairs {+z_i,-z_i}
of the verified k-dimensional direction code. The assignment is explicit: choose
one representative of each antipodal coefficient tuple by first-nonzero sign,
sort these tuples lexicographically, and match them to class order. For every
selected signed unit mother head u in class i, remove u from the equator and
insert the two unit points

```text
(sqrt(3)/2 * u, +z_i/2), (sqrt(3)/2 * u, -z_i/2).
```

Retain every other mother point at the equator, and add the separately verified
unit pole code. For A2, E6 and E7 tails, the exact displayed linear equations give
subspaces of dimensions two, six and seven; use any orthonormal identification
with those Euclidean spaces. Thus the ambient dimension is at most 48+k.

All inner products are bounded by 1/2:

1. In one class, different non-antipodal heads give at most
   (3/4)(1/3)+1/4 = 1/2. Opposite heads are also safe.
2. The two copies of one signed head have product 3/4-1/4 = 1/2.
3. Across classes, distinct mother heads have product at most 1/2 by the P48
   minimum. Distinct direction pairs give at most 1/2, so the mixed product is
   at most 3/8+1/8 = 1/2.
4. Retained equator/mixed products are at most sqrt(3)/4, strictly below 1/2.
5. Mixed/pole products are at most 1/2 by Cauchy–Schwarz. Equator/pole pairs are
   orthogonal; equator/equator and pole/pole pairs satisfy their own kissing bounds.

The direction and pole codes need not be mutually compatible. In 60D they have
840 and 841 points respectively. With L selected mother lines and P poles,

```text
N = 52,416,000 - 2L + 4L + P = 52,416,000 + 2L + P.
```

The 60D witness has L = 3,105,435 and P = 841. The entire 58-million-point code is
specified by this recipe; it is not expanded and compared pairwise.

## Discovery and reproduction

The initial isometries are the published image selections. Starting with their
full union, assign each line to its first selected image. A candidate may be
inserted when all its conflicting retained heads are absent or can be moved to
another image that contains the same line and remains compatible with its added
heads. Reassignment preserves the occupied mother lines; insertion enlarges the
union. Every final destination class is checked from its actual contents.

A deterministic pool comes from 3,500 products of 30 published generators,
using NumPy `default_rng(20261108)`. It retains 22,643 novel candidates with at
most three conflicts against the complete base. Candidate and class order,
three-pass limit and tie-breaking are recorded in the supplied scripts. This is
a restricted heuristic search, not an exhaustive shell, class or ownership search.

To recreate the search in a separate output directory:

```sh
python discovery/sieve.py /path/to/output
python discovery/augment.py /path/to/output
```

The second command uses only the initial transforms and IDs from `heads*.npz`,
not their final masks or extra points, and rebuilds ownership from scratch.
To skip pool generation, copy `discovery/low_conflict_pool.npz` into the output
directory before running `augment.py`. Compare arrays rather than ZIP metadata.
Search history is not trusted by either final geometry verifier.

## Provenance and scope

Inputs are pinned to repository commit
[52eac7a24a5912ac5262733165f55abfe6869aa8](https://github.com/alexlegeartis/KissingNumbers/commit/52eac7a24a5912ac5262733165f55abfe6869aa8).
`source/SOURCE_PROVENANCE.json` records head inputs and initial image selections;
`source/TAIL_SOURCE_PROVENANCE.json` records direction and pole inputs, source
formats and the MIT license. The 840-point construction is from
[Takhanov](https://arxiv.org/abs/2609.09179). All prior witness files and verification programs are unchanged. Historical
README and index text, with their document checksums, identify supersession.
The manuscript is unchanged.

Larger mother classes, broader candidate pools and different tail codes may
improve these constructions further. No source-pool optimum, unrestricted
optimality or worldwide priority is asserted.
