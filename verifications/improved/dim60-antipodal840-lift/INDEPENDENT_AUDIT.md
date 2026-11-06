# Independent audit of the 60-dimensional lift

Audit date: 6 November 2026. Verdict: the proposed geometry and count are valid conditional on the three input certificates: the P48 mother shell and disjoint classes, the antipodal 840-point direction code, and the independent 841-point pole code. This audit checks the lift proof and stored class counts independently; it does not claim a fresh replay of all three inputs.

## Exact geometry

Normalize the P48 mother shell C to unit norm. For each selected disjoint class of lines, let B_i contain both signs of its heads. Distinct heads in a single B_i have inner product at most 1/3. Distinct heads anywhere in C have inner product at most 1/2. Choose one unit representative u_i from each of 420 distinct antipodal pairs of the direction code, so |u_i dot u_j| <= 1/2 for i != j. Let Y be an independently chosen 841-point unit kissing code in R12.

Retain (x,0) for heads outside the selected classes. Replace each selected signed head x in B_i by the two points (sqrt(3)/2 x, +u_i/2) and (sqrt(3)/2 x, -u_i/2). Append (0,y) for every y in Y.

All points have unit norm. Pair bounds are:

- Two retained equatorial points: at most 1/2.
- Equatorial versus mixed: at most sqrt(3)/4 < 1/2, because the selected head itself is no longer on the equator.
- Mixed points in different classes: at most 3/8 + 1/8 = 1/2. Disjointness of the classes is essential here.
- Mixed points in one class with distinct signed heads: at most 1/4 + 1/4 = 1/2.
- Two mixed points with the same signed head: the tail signs differ, giving 3/4 - 1/4 = 1/2.
- Mixed versus pole: absolute inner product at most 1/2 by Cauchy-Schwarz alone.
- Two poles: at most 1/2 by the pole-code certificate.

Thus Y need not coincide with the direction code, contain it, or be antipodal. The 841 poles and 840 directions are independent ingredients. No new mixed-pole pair computation is required beyond their unit norms. Distinct supports/norm splits and distinct selected heads prevent duplicate points.

## Independently checked arithmetic

The stored family_refilled.npz has 577 nonincreasing class sizes. The sums of its largest 378 and 420 entries are respectively 2,775,751 and 3,062,990. Every stored class size agrees with the number of retained base indices obtained by unpacking retained_bits in the default big-endian bit order. The file SHA-256 is d2f6afe15f1ad05b08832e4e9dcbf0a7a1159e232ab63ad30259b98ce60e0396, matching the supplied verification receipt.

For S=3,062,990 selected lines there are 2S removed equatorial points and 4S added mixed points. The resulting count is

52,416,000 - 2S + 4S + 841 = 58,542,821.

Against the current repository's 57,968,343 this is an increase of 574,478, exactly twice the 287,239 additional retained lines. There is no +84 pole gain against this baseline, because it already uses 841 poles.

Selecting whole already-disjoint classes from the 577-class certificate is safe. Merely summing untrimmed image sizes would not be safe. This audit has checked the stored count/bitmap agreement; the construction worker is responsible for replaying head-class geometry and disjointness.

## Sources and comparison scope

- https://arxiv.org/html/2609.09179v1 : Theorem IV.1 explicitly states an antipodal 840-point code; Appendix Lemma .6 proves antipodality using all sign lifts. This is the correct source, not an assumption that every 840-point code is antipodal.
- https://github.com/alexlegeartis/KissingNumbers/blob/main/RESULTS.md : fetched directly through the GitHub connector on the audit date; 60D is 57,968,343 and still uses 378 direction pairs.
- https://github.com/alexlegeartis/KissingNumbers/blob/main/verifications/improved/dim49-63-p48-caps/README.md : fetched directly; confirms the 7,783-line base, 577 classes, 756-point direction code, independent 841-point pole code, and literature dependence of the mother-shell count/minimum.
- https://arxiv.org/html/2611.01906v1 : November 1 paper reports the weaker 57,573,374 using 378 pairs and 756 poles. Its repository README also still reports this value when inspected.

A current web search did not reveal a stronger public 60D construction or this 840-direction substitution. That is a bounded source check, not proof of worldwide novelty. The earliest +84 observation against the November 1 paper is superseded by the stronger repository baseline and must not be presented as the final improvement.

## Subsequent joint-coverage completion

The updated 58,543,525-point result was independently audited after the initial publication. Fresh full-union counts replayed all 37 recorded exchanges: 3,062,990 lines became 3,063,342. The packed retained memberships cover the final selected union exactly once. A separate endpoint swap calculation using set intersections, rather than the optimizer's owner-vector shortcut, found maximum one-swap gain zero. This is a weak local optimum: zero-gain moves remain, and no global optimality follows. The full geometric verifier was rerun successfully, including 420 isometries, all 30,283,653 base pairs and canonical 48-coordinate disjointness. The optimized420.npz SHA-256 was b3f9c757533949c90601ae0de2adb4386bdf20371da1b852740595bc1da5cdbf. The inherited lattice and independent pole premises above are unchanged.
