# Nine-coordinate affine caps in 61–66 dimensions

This is an **append-only** improvement to six certified kissing-number lower bounds. No old points are removed. Count the finite new layer independently with `OPENBLAS_NUM_THREADS=1 python verifications/improved/dim61-66-nonary-affine-caps/verify.py` from the repository root. The discovery search is not needed to replay the certificate.

| Dimension | Pinned incumbent | New lower bound | Gain |
|---:|---:|---:|---:|
| 61 | 88,979,562 | **89,072,802** | 93,240 |
| 62 | 122,543,436 | **122,659,076** | 115,640 |
| 63 | 184,415,700 | **184,417,784** | 2,084 |
| 64 | 335,749,696 | **335,789,656** | 39,960 |
| 65 | 335,817,922 | **335,857,882** | 39,960 |
| 66 | 342,291,652 | **342,297,172** | 5,520 |

## Construction

In the incumbent's scaling, squared norm is 64 and the dot-product limit is 32. Each new point has nine occupied coordinates of equal amplitude 8/3. On each support use the 40 binary sign words obtained by puncturing the Gray image of the quaternary pentacode ([Conway–Sloane, 1994](https://neilsloane.com/doc/Me185.pdf)). The verifier directly enumerates the cyclic shifts and negations of the four quaternary seeds, applies the Gray map, punctures, and checks all 780 word pairs and minimum distance three. The sign code is a finite certificate, not a bound assumed from the literature. The squared norm is `9*(8/3)^2=64`; within one support the dot product is at most `(9−6)*(8/3)^2<32`. Across different new supports the checked intersection at most four gives at most `4*(8/3)^2<32`.

The inherited weight-sixteen points use amplitude two. The verifier reconstructs every support, independently of signs, and checks each new/old intersection at most six. This gives `6*2*(8/3)=32`. For the amplitude-four layer we check no new support contains one of its four-coordinate blocks. In 61–65D an affine SQS has blocks with zero coordinate XOR, so pairwise distinct XORs among the 36 pairs in a new support certify this condition. In 63D the test uses the incumbent's explicit permutation of the affine SQS. In 66D the verifier regenerates the pinned Netto four-support stock and compares its blocks directly. Thus the new/weight-four dot product is at most `3*4*(8/3)=32`. Any axis has dot product at most `8*(8/3)<32`. The dense layers have amplitude one except for the small tail coordinates in 61D and 62D and the radical body scaling for one 61D latitude. On nine occupied coordinates the 62D maximum absolute coordinate sum is at most 11; in 61D the two latitudes give at most 10 and `9*sqrt(17)/4`. All are less than 12, which suffices.

For **63D**, take the specified signs with zero signed sum on the intersection with the exact pinned section octad. Per-support permutation of the same 40 signs gives the stored number of balanced choices. The ambient 64-dimensional inequalities imply the same inequalities in the section. For **65D**, use the 64D new points at zero tail. The incumbent's owner lifts are over the dense and four-support families already covered by the bounds above, with ordinary factor sqrt(3)/2. Both new axial poles are orthogonal to this layer.

## Inherited proofs

- [61/62/66 affine sign reconstruction](../dim61-62-66-affine-reconstruction/README.md), with the [exact dense, tail and sparse layers](../dim61-67-compatible-sign-families/README.md).
- [63/64/65 affine sign reconstruction](../dim63-65-affine-reconstruction/README.md), including the frozen 63D octad and SQS relabeling. The [65D owner lift](../dim64-65-signed-parity/construction.py) is unchanged.

The verifier pins the incumbent support archives, certificates, coordinate conventions and the Netto constructor by hash. The earlier configurations' own geometric validity and counts are inherited, not replaced by this checker. No optimality claim is made for the search, the new layer, or these dimensions.
