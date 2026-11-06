# An 840-direction cap lift in dimension 60

**Exact replay gives K(60) >= 58,543,525**, an increase of **575,182** over the repository's pinned 57,968,343 comparison. This combines existing constructions; it does not claim a new general lifting theorem or establish worldwide priority.

The missed ingredient was the distinction between the largest known kissing code and a large antipodal code. The 841-point 12D record is asymmetric, but Takhanov's separate 840-point quaternionic construction is antipodal. It supplies 420 direction pairs, instead of the 378 pairs from the 756-point Coxeter–Todd shell used by the existing cap package. The 841-point code remains useful independently as pure poles.

## The construction

Use the published P48p minimal shell of 52,416,000 norm-six vectors. Select 420 complete class images by the saved maximum-coverage exchanges, then assign every union line to one containing image. They contain **3,063,342 unoriented lines**. Each class has absolute normalized products at most 1/3; different mother lines have absolute normalized products at most 1/2.

Enumerate the 420 antipodal direction pairs {±z_i} in the exact 840-point tail. Match the sorted pairs to the selected classes in their saved order in optimized420.npz. For each selected mother line {±u}, with u normalized, replace its two equatorial points by the four points

    (epsilon * sqrt(3)/2 * u, eta/2 * z_i),  epsilon, eta in {−1,+1}.

Retain all other equatorial minimal vectors, and add the separate 841-point unit code as pure poles (0,y). Therefore

    52,416,000 − 2S + 4S + 841
    = 52,416,000 + 2 × 3,063,342 + 841
    = 58,543,525.

The pair-type proof is short. Within a class, the worst product is (3/4)(1/3)+1/4=1/2. Across classes it is (3/4)(1/2)+(1/4)(1/2)=1/2. Distinct points on one selected mother line have maximum product 3/4−1/4=1/2. Equator–cap products are at most sqrt(3)/4. Cap–pole products are at most 1/2 by Cauchy–Schwarz, regardless of the relation between the direction code and the pole code. Thus the two tail codes need not coincide or be mutually compatible.

## Exact finite evidence

Run:

    OPENBLAS_NUM_THREADS=1 python verifications/improved/dim60-antipodal840-lift/verify.py

The replay uses Python and NumPy, without an optimizer. It checks:

- All 352,380 pairs in the reconstructed 840-point direction code in Q(sqrt(2)), and all 420 antipodal pairs. Coordinates are stored as (A+B sqrt(2))/4.
- All 353,220 pairs of the 841 pure poles, treating the published decimal coordinates as exact rational rays and using Python arbitrary-precision integers.
- All 30,283,653 pairs of the 7,783-line base class, using bounded exact int64 products with the published Gram matrix.
- All 420 selected isometries and retained vector norms.
- Disjointness of all 3,063,342 mother lines by comparing their complete canonical 48-coordinate representations, not probabilistic hashes.
- The algebraic bounds for every type of pair in the lifted construction.

The mother's minimum six and shell population 52,416,000 remain explicit published lattice inputs. The 58-million-point construction is specified by a finite recipe, not expanded and compared pairwise. The certificate must not be described as a fresh enumeration of that entire shell or all its point pairs.

`verify_heads.py` without --full also verifies a historical-file-independent fallback with 3,035,307 lines and 58,487,455 points, using only generator-word transforms. The stronger certificate uses the historical transforms too, each independently checked to preserve the Gram matrix.

## Sources and comparison

The repository inputs are pinned to commit **11dc33180eb63bcf255f81a9a4cc95f187bdd174**. The historical binary was retrieved by filtered git clone after the contents connector returned an empty payload; its Git blob hash is 11cb2a525c813b2d605672d479c5d6c5df3f12d5. SHA-256 hashes and file sizes are in SOURCE_MANIFEST.json.

- P48 cap source: https://github.com/alexlegeartis/KissingNumbers/tree/11dc33180eb63bcf255f81a9a4cc95f187bdd174/verifications/improved/dim49-63-p48-caps
- Pinned comparison registry: https://github.com/alexlegeartis/KissingNumbers/blob/11dc33180eb63bcf255f81a9a4cc95f187bdd174/RESULTS.md
- Takhanov's 840-point quaternionic construction: https://arxiv.org/abs/2609.09179
- Antipodal block lift and independent November1 comparison: https://arxiv.org/html/2611.01906v1
- Earlier general lift: https://arxiv.org/html/2609.35051v1
- 841-point pole coordinates: the exactly verified source file preserved in this package; source provenance is shared with the repository's 12D paper audit.

The November1 paper's 60D value is 57,573,374, weaker than the repository comparison. Both use 378 direction pairs. The present improvement is over both checked sources, not a comprehensive proof that no other unpublished or unindexed construction is stronger.

## Transferable lesson

For a mother configuration with suitable antipodal block classes, select the direction code by its number of antipodal pairs and select the pole code independently by total cardinality. A new low-dimensional direction code can improve a much larger construction without improving that low dimension's unrestricted record. Here the known 840-point code, although smaller than 841, activates 42 more head classes and changes the 60D bound substantially.

## Maximum-coverage completion

The initial published composition (commit f0083954ce63e292b0d6024d73e35cb2486716e2) used 3,062,990 lines and gave 58,542,821 points. Refilling its same 420 images gains nothing. However, ranking pre-trimmed classes misses overlaps assigned to omitted classes. Thirty-seven deterministic one-class exchanges among the 577 complete image sets increase the union by 352 lines, hence 704 points. Each line is assigned to the first selected image containing it, preserving the within-class bound. The final masks and isometries are in optimized420.npz; coverage_exchanges.json records the exchanges. verify_coverage.py independently checks the selected geometry and disjointness without trusting the optimizer. This is a one-exchange local optimum in the supplied finite image family, not a global optimum. The earlier head certificates remain as historical evidence.
