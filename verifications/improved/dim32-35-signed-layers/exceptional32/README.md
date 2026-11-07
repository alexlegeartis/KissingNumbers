# A 352,832-point kissing configuration in dimension 32

Checkpoint: **7 November 2026**. The frozen integer certificate proves
**tau(32) >= 352,832**, an increase of 32 over the reviewed 352,800-point
repository construction and 512 over the pinned 352,320-point Lin result.
Every point has squared norm 32, and every distinct pair has inner product
at most 16. Scaling by 1/sqrt(32) gives maximum inner product 1/2. The roots
span all 32 coordinates. No global optimality or worldwide-priority claim is made.

## Reproduce without a solver or downloads

Python 3.10+ and NumPy suffice. From this directory, using an available CPU
and keeping at least one logical CPU free:

```sh
taskset -c 0 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python replay.py
```

The parent `../verify.py` checks the complete package manifest, copies this
folder to a temporary directory and runs the same replay. Thus the saved
receipts are evidence, never inputs to the mathematical conclusion.
`replay.py` fails under Python `-O`, which would remove assertions.

- `verify.py`: original mask/code proof, applied to the new certificate.
- `verify_coordinates.py`: original independent coordinate proof, expanding
  every dense, middle and root point and checking exceptional interactions.
- `independent_audit.py`: a fresh independent implementation, with no imports
  from either preceding checker or the search. Materializes all 352,832
  vectors and checks their norms and uniqueness. Only its default paths and
  usage line were adapted for this self-contained package.
- `verify_delta.py`: verifies the complete 380-point removal/412-point addition
  listed in `delta.json`, plus unchanged dense generators and full supports.
- `verify_pool_bound.py`: constructs exact conflict-clique covers proving the
  1,664-point capacity for these 25 exceptional supports only.
- `negative_controls.py`: four corruptions rejected by each original checker.
  The independent audit rejects nine different data/geometric corruptions.
  `replay.py` additionally verifies rejection of an incorrect expected hash.

`verification.json`, `coordinate_verification.json`, `report.json`,
`pool_bound.json` and `negative_controls.json` are replay receipts.
`evidence/report.json` and `evidence/hash_negative_control.log` preserve the
pre-publication independent audit and hash-rejection evidence.

Candidate SHA-256:
`69c3f368a4aae5956ca75a8569f761224148831e8a85edd5972a0b3b057af646`.

Retained control SHA-256:
`3a96dc79d21d07557a221f2cdd6664d2d5270e93d7c98711d7b962a2ba2b41d2`.

`control32_352800.json` is an exact byte copy of the unchanged historical
`../source/classic/certificate32.json`. The candidate's inherited `comparison`
field is historical metadata, not a fresh literature survey.

## Complete pair coverage

The four families contain 131,072 dense vectors, 218,112 full-octad vectors
(1,704 supports times 128 even sign patterns), 1,664 exceptional octad vectors
and 1,984 roots. All coordinates are integers.

- Dense vectors are signs of magnitude 1 from the supplied binary [32,17,8]
  generator. Enumerating all 131,072 different words proves rank 17 and
  minimum nonzero weight 8. Linearity bounds every dense/dense product by 16.
- Full octads have eight coordinates of magnitude 2. All 8,128 unordered
  pairs in the 128-pattern even sign cube have product at most 16. The
  1,450,956 distinct support pairs have intersection at most four, bounding
  every cross-support product by 16.
- The 832 exceptional lines expand into 1,664 distinct vectors on 25
  supports. The independent audit directly checks all 1,383,616 unordered
  exceptional products, including same-support and antipodal pairs.
- All 42,600 full/exceptional support pairs intersect in at most four places,
  bounding every mixed product by 16.
- Roots have two coordinates with independently chosen signs of 4. All
  1,967,136 unordered root pairs are checked directly, with maximum 16.
- The remaining mixed products have pointwise shape bounds: dense/octad
  at most 8*1*2 = 16, dense/root at most 2*1*4 = 8, and octad/root at most
  2*2*4 = 16.

This covers all ten unordered family-pair types. It is a complete finite
proof, not a claim to multiply the entire 352,832-square Gram matrix.
All products use int16 or wider; intermediate magnitudes are at most 32.
The independent audit reconstructs the old control and applies the same
proof to it. Both full coordinate arrays have the claimed number of distinct
rows and squared norm 32 throughout.

## Exact change and restricted capacity

Only nine exceptional supports change. Their component grows from 544 to
576 points. The other support-graph components keep 192, 384, 384 and 128
points, giving 1,664 exceptional points in total. Dense vectors, full octads
and roots are unchanged. The search allowed individual signed points, without
imposing antipodality, affine fibers or a parity restriction. The solution
happens to be antipodal and can therefore use the historical line format.

The 25 supports partition into eight triangles and one isolated support.
For octad supports meeting in five, six or seven coordinates, fix the signs
on their intersection, match the remaining sign choices, and flip one
coordinate outside the intersection on either side. This partitions all
512 signed vectors on the pair into 128 conflict cliques of size four.
Every within-clique product is greater than 16, so at most 128 vectors can
be selected from either pair of supports. Summing the three pair bounds
in a triangle gives 192. A single octad cube has capacity 128 by its
sign-flip perfect matching. Hence 8*192 + 128 = 1,664, attained here.
`verify_pool_bound.py` constructs all 24 pair covers, verifies disjointness,
complete coverage and every conflict in integer coordinates.

This upper bound concerns only these same 25 supports with the other layers
fixed. It is not an upper bound for other supports, other shapes, or tau(32).
An exploratory search on 171 compatible octad supports found no additional
gain; that experiment is not an optimality proof or a premise of this package.

## Comparisons and attribution

The 7 November read-only review resolved source branches and read the files at
exact commits. `evidence/public_comparison.json` records URLs, timestamps and
excerpts. Corresponding source snapshots are retained as `.txt` or JSON.

- [Repository RESULTS at b9164ae](https://github.com/alexlegeartis/KissingNumbers/blob/b9164aebd6b0e09105b398ed7d9a593bf52a2df5/RESULTS.md)
  lists 352,800. Its certificate is structurally identical to the control;
  the saved tool snapshot normalizes line endings but the control retains
  its original bytes.
- [Lin README at f367a9c](https://github.com/QH-Lin/kissing-number-32/blob/f367a9c5a051c0d6821ce71324b9dfa1ead76c2e/README.md)
  gives 352,320. Its README, rather than the stale repository description,
  is the comparison source.
- [Cohn's table](https://cohn.mit.edu/kissing-numbers/) could not be freshly
  retrieved (timeout/HTTP 502). No fresh Cohn value is asserted.

Thus the verified gains are +32 over the reviewed repository and +512 over
Lin. These source checks do not certify worldwide priority. The repository's
SOTA label is relative to its named comparisons and this stated review scope.

The public support codes are due to Qihong Lin. The dense generator is
Cheng-Sloane's [32,17,8] code, with provenance retained in `../source/public/`.
The underlying layers follow [Edel-Rains-Sloane](https://arxiv.org/abs/math/0207291);
signed replacements also appear in [Yang et al.](https://arxiv.org/html/2609.35051v1).
All historical certificates and the manuscript remain unchanged.
