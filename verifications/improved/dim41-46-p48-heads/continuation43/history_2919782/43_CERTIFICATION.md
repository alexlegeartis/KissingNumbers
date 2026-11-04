# Dimension 43: current certificate and historical gap

**28 October 2026:** the registry now contains **2,919,782** points.
The continuation deletes 46 points from the preceding 3,408-point tiny layer
and inserts 85, for a gain of **39**. The current 3,447 explicit tiny points
are in `partial43/explicit_points_43.npz`; the predecessor and exact delta
are retained in [continuation43/](continuation43/README.md).

The 118 complete families and five pure heads retain their previous inputs.
Both `partial43/verify.py` and `continuation43/verify_coordinates.py` pass
6,370,783 angle conditions and independently replay 22,001 short-difference
refinements. The second checker uses the rank-40 coefficient metric and
integer square-root cutoffs. Complete-family populations still inherit the
modular identities and published P48n minimum/design inputs; not every large
shell is freshly enumerated. The exact coordinate delta needs no optimizer.

The historical missing 2,911,990-point target remains unnecessary for this
construction and its 2,919,743-point predecessor.

## Earlier checkpoints

- 2,878,155: recovered seed in `data/old43_seed.npz`. `certify43.py` freshly
  enumerates its 13 large shells and checks its angles. The recorded report
  is `certification/recovered43.json`.
- 2,911,990: historical writeup target. The supplied writeup referred to
  `kissing41-46_exact_certificates.zip` and `kissing43_extension_2911990.zip`,
  which were unavailable at the first integration. That target was excluded.
- 2,918,396 and 2,919,686: later complete-family constructions in the source
  research. Partial-family trades then gained 57 points, producing the
  then-current 2,919,743-point certificate, now retained in the continuation package.

The historical missing-certificate finding remains about the older target.
It is no longer an outstanding input for the registered 43D construction.
