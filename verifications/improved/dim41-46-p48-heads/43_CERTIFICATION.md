# Dimension 43: current certificate and historical gap

**24 October 2026:** the registry now contains **2,919,743** points, certified
by `partial43/verify.py`. This supplied construction supersedes the earlier
uncertified 2,911,990-point writeup target. It does not use the missing target
archives or the seven extension points named in that writeup.

The current construction has 118 complete body families, 3,408 explicitly
selected tiny-body/head pairs and five pure heads. Fresh integer replay
passed 6,232,333 angle checks and 21,825 short-difference refinements. Exact
parent embedding, the full rank-40 complement, gluing, shortest bodies and
the rank-three head projector replay. Complete-family populations inherit
the modular identities and published parent assumptions described in the
[package README](README.md); not every large shell was freshly enumerated.

## Earlier checkpoints

- 2,878,155: recovered seed in `data/old43_seed.npz`. `certify43.py` freshly
  enumerates its 13 large shells and checks its angles. The recorded report
  is `certification/recovered43.json`.
- 2,911,990: historical writeup target. The supplied writeup referred to
  `kissing41-46_exact_certificates.zip` and `kissing43_extension_2911990.zip`,
  which were unavailable at the first integration. That target was excluded.
- 2,918,396 and 2,919,686: later complete-family constructions in the source
  research. Partial-family trades then gained 57 points, producing the
  current supplied 2,919,743-point certificate.

The historical missing-certificate finding remains about the older target.
It is no longer an outstanding input for the registered 43D construction.
