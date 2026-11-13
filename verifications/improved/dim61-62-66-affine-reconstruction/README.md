# Affine sign-subcode extensions in dimensions 61, 62 and 66

Reconstructed on 13 November 2026 from the published compatible-sign families.
These are explicit, certified lower-bound constructions, not claims of optimality.

| Dimension | Previous repository construction | This construction | Added points | Certificate and proof checks |
|---|---:|---:|---:|---|
| 61 | 88,975,146 | **88,979,562** | 4,416 | [certificate](dim61-62/certificate.json), [independent replay](dim61-62/independent_verification.json) |
| 62 | 122,539,020 | **122,543,436** | 4,416 | [certificate](dim61-62/certificate.json), [independent replay](dim61-62/independent_verification.json) |
| 66 | 342,271,364 | **342,291,652** | 20,288 | [certificate](dim66/certificate.json), [independent replay](dim66/independent_verification.json) |

## Idea: keep partial sign classes instead of rejecting entire supports

The sparse layer has weight-16 supports with entries ±2, hence squared norm 64.
An ordinary extended-Hamming sign class has 2,048 points. A candidate support
that overlaps an existing support in ten coordinates need not be discarded:
retain an affine subcode of signs satisfying extra parity equations. On each
critical ten-coordinate intersection we impose a common weight-eight parity
check with opposite phases. Thus at least one sign disagreement occurs on the
intersection and the inner product is at most 4(10−2)=32. Intersections of
size at most eight already satisfy the same bound regardless of signs.
The linear sign codes have minimum Hamming distance at least four, giving
inner product at most 64−8·4=32 within each class.

To accommodate these checks, a few old Hamming classes are replaced by other
full 2,048-point Hamming classes on the same supports. This changes signs but
loses no points. There are 25 new subcodes and one modified full class for
61/62D, and 28 new subcodes and 104 modified full classes for 66D.

In 61/62D the supports avoid coordinates 24, 42, 53, 55. Odd sign parity and
orthogonality to the dense binary-code generators preserve compatibility
with the dense layer. Each additional sparse point transfers to both dimensions
using the original dimension-transfer construction. In 66D all parent supports,
including the 3,200 supports with nonzero tail coordinates, are included in
the checks; the independent verifier reconstructs these from the cyclic73 archive.
All unchanged layers and transfer geometry are those of the
[original package](../dim61-67-compatible-sign-families/README.md).

## Reproducibility and scope

From the repository root (Python 3 and NumPy):

```sh
python verifications/improved/dim61-62-66-affine-reconstruction/verify.py
python verifications/improved/dim61-62-66-affine-reconstruction/verify.py --independent
```

Both entry points also work from any working directory. The fast verifier checks
exact parity certificates and every affected support pair. The independent
verifier instead filters all 65,536 possible sign words and directly computes
critical restricted-sign inner products, without importing discovery helpers.
It also confirms that the embedded baseline family equals the published file.
Assertions are essential; optimized Python (`-O`) is explicitly rejected.

The baseline replay reports record freshly completed checks of the original
construction, including the scalar 2^28-word binary-code replay. The 66D
baseline replay includes 578,935,378 parent-support pair checks and the Netto
layer checks. Those baseline routines are the published implementations,
not independently rewritten geometry proofs. To rerun all baseline checks:

```sh
python verifications/improved/dim61-67-compatible-sign-families/verify.py
```

The certificates specify each support as an integer bitmask, local parity-check
rows, an affine sign shift, and every added class size. Bit order is increasing
support-coordinate order. `context.json` retains the corresponding discovery
context; the independent verifier does not trust its 66D parent list without
reconstructing and comparing it. Source baseline commit:
`69649636c78c8166fe277d7118dc21435d6765d0`.

`final.py` derives registry counts from the certificates. The independent entry points are `dim61-62/verify_independent.py` and `dim66/verify_independent.py`.
