# Exact 39-dimensional lower bound: 769,852

The producer's complete integer verification and the independent direct-integer audit both pass. Four arbitrary-precision witness replays and eight rejection controls also pass. This is a constructive lower bound, not a claim that the kissing number has been determined exactly or that worldwide priority has been exhaustively checked.

## Construction

The sparse layer contains 441,244 vectors of squared norm 32:

- 3,424 weight-eight supports, carrying all 128 even-parity sign patterns with nonzero coordinates of magnitude 2: 438,272 vectors.
- The retained signed roots, with two coordinates of magnitude 4: 2,908 vectors.
- Four heavy stars, each contributing 16 vectors: 64 vectors.

This layer is unchanged from the preceding certificate. Its support replacements, omitted root pairs, heavy stars and original support data are preserved in `source/`.

The dense layer has 10,269 integer representatives and a 32-element coordinate-sign group, yielding 328,608 points. Each integer direction is divided by its own exact Euclidean norm. Thus the total is 441,244 + 328,608 = 769,852, an increase of 32 over the preceding 769,820-point construction.

## Exact proof

Let L = 99,999,996 and U = 100,000,004. Every dense integer direction has norm in [L,U]. The verifier checks the sign group, its complete action and every relative dense product, excluding only the identity action on a point with itself.

There are 1,687,391,811 checked relative products. The largest is 4,999,991,179,259,176, below floor(L²/2) = 4,999,999,600,000,008. Consequently every normalized dense pair has inner product at most 1/2.

The verifier also checks 35,161,056 actual middle-support bounds, using the group projection parity to cover precisely the required sign possibilities. The minimum exact squared middle slack is 23,357,272,208. All root and heavy bounds pass with strict positive integer slack. The sparse verifier checks the support packing, retained roots, heavy stars, all sparse pair classes and exact counts. Together these checks cover every pair in the full configuration. Arithmetic bounds are checked before using signed 64-bit arrays; no floating-point tolerance or optimization assumption is used for acceptance.

Producer receipt: `orbit_verification.json`.

Frozen certificate SHA256: `dee846d08a9516ca12a55fcf89bde154e529ba6302d7146487bed92d18a58300`.

## Discovery and scope

A numerical 32-extra-orbit pose was pruned to 29 extra orbits by exhaustively pricing all 4,960 deletion triples. All dense representatives were then jointly deformed. The preceding dense configuration is not assumed to remain fixed. Conservative angular screening and a checked extrapolation only accelerated the numerical search; the final exact proof enumerates every required pair independently of that screening.

The separate five-support trade search found no sparse support-count improvement when deleting at most five supports. That local negative result is not used in this proof.

This package publishes the frozen certificate and preserves every previously published dense certificate.

## Independent coverage

The separately implemented direct-integer route checks all 1,687,391,811 dense relative products without a transform, all 8,859,271,680 permitted middle sign sums, 21,030,912 dense/heavy checks, all retained root coordinate-pair bounds, and every sparse pair class. It reproduces the exact dense maximum and its witness. Its coordinate guard is derived from U; all signed 64-bit accumulation and squaring bounds are checked. The untouched independent audit archive is `independent_original.tar.gz`, with SHA256 `e1f83f77997aef9e6613fb0d9db4a2a70de090039b0a249f4d5f258f306e5431`. This is a finite computational certificate, not proof-assistant formalization.
