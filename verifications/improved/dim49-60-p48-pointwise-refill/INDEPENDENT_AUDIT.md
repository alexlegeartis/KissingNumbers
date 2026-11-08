# Independent finite audit

All twelve final head witnesses passed the separate global-coordinate checker.
Its per-witness receipts pin the exact certificate, source and verifier bytes.
The final 60D witness is `cc3a2e07fbe3f5f91c33d61ba69102057c5fd8f90dc509d3c5d83d29adac9b25` and contains 3,105,435 distinct
antipodal mother lines, giving 58,627,711 points after the verified lift.

The independent implementation does not reuse the discovery conflict lists,
owner dictionary, reported gain or asserted objective. It reconstructs each
final class, checks every isometry and norm, checks all added-head pairs against
its actual final class, and compares full canonical global coordinates for
ownership. Across the twelve witnesses it checks 655,140,587 added/retained
or added/added pairs, in addition to the 30,283,653 base pairs.

Its binary64 arithmetic is exact: each input is an integer and each product and
absolute partial-sum bound is strictly below 2^53. The separate producer checker
uses bounded int64 products in local coordinates. Their distinct coordinate
frames and arithmetic implementations guard against a shared cached-conflict or
ownership interpretation error.

The tail checker verifies all 751,195 direction pairs and 353,220 independent
pole pairs. Dimension reductions for A2/E6/E7 use exact linear relations. The
841 printed pole coordinates are treated as exact rational rays, each normalized
by its own norm; no numerical tolerance or rounded unit-norm assumption is used.

Both head and tail proofs retain the published identification, minimum six and
population 52,416,000 of the P48p Gram lattice. These are not fresh mother-shell
or full lifted-code enumerations. The twelve mutation controls are listed in
`audits/negative_controls.json`; run `test_certificate.py` to repeat them.
