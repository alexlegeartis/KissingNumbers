# Signed-owner exchange transported across P48p tail classes

The jointly improved signed stock in 49D gives a shared pool of 343 literal norm-six donors. Transport them into the existing tail classes using the existing P48p isometries. For each class optimize the subset of these donors against its actual existing signed owners: one point gained for each inserted donor, one point lost for each incumbent owner that blocks any chosen donor. This overlap-aware optimization is a maximum-closure/min-cut problem because the common donors are jointly compatible. Global ownership exclusions prevent assigning a mother direction twice.

| Dimension | Prior lower bound | New lower bound | Gain |
|---:|---:|---:|---:|
| 50 | 52,462,872 | **52,462,974** | +102 |
| 51 | 52,509,740 | **52,509,940** | +200 |
| 52 | 52,603,402 | **52,603,800** | +398 |
| 53 | 52,728,054 | **52,728,705** | +651 |
| 54 | 52,976,623 | **52,977,769** | +1,146 |
| 55 | 53,393,546 | **53,395,469** | +1,923 |
| 56 | 54,263,398 | **54,266,509** | +3,111 |
| 57 | 54,760,693 | **54,764,394** | +3,701 |
| 58 | 56,270,076 | **56,274,806** | +4,730 |
| 59 | 56,952,414 | **56,957,415** | +5,001 |

## Exact replay

Run `python verify.py` (NumPy). If the public parent certificate directories are not siblings here, supply `--pointwise` for the existing dim49-60-p48-pointwise-refill and `--signed` for dim49-60-p48-signed-refill. Checksums pin those two parent class-geometry inputs and the common donor pool. The compact gzip witnesses specify only which signed donors each class receives. The verifier reconstructs the previously certified class ownership, derives every blocker afresh by exact products, checks donor norms, new/new compatibility, class isometries, and complete-coordinate global ownership. No optimizer, confusion graph, or geometric tolerance appears in the replay. All changed-class products use bounded integer-exact binary64 accumulations. The unchanged old components, lattice identification, equatorial kissing count, and tail/pole geometry come from the pinned parent certificates.

For each selected unit mother point u in a particular tail direction line z, replace its old equatorial point by `(sqrt(3)u/2, +z/2)` and `(sqrt(3)u/2, -z/2)`. Within a class the old and new distinct owners have normalized product at most 1/3, so same-tail products are at most 1/2. The two copies of the same owner have product 1/2. Different owners in different classes have product at most 1/2 in the body, and distinct tail lines at most 1/2 in absolute value in the head, also giving ≤1/2. Entirely distinct signed ownership is necessary across classes and is rechecked. All unselected mother points and the prior pure pole code remain unchanged. Counts follow as 52,416,000 plus selected signed owners plus the original number of poles.
