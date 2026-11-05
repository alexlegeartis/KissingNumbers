# Individually selectable rational group images

`search.rational_image_exchanges` takes a verified integer spherical code with
positive diagonal metric and exact rational isometries. A generator is
`{"numerator": [[...]], "denominator": d}`. Every generator is independently
checked by the exact identity A^T W A = d^2 W. `words` is a finite list of
lists of generator indices; a word `[0, 1]` applies generator 0 and then 1.
The empty word is the identity. If omitted, each generator is used once.
Noncommuting words can therefore be searched without declaring a group finite.

All words are represented over one exact common denominator. The original
configuration and the images are deduplicated and candidates are individually
selectable: an image need not be taken as a complete orbit. Incumbent points
are also independently removable. Shared source deletions are charged only
once. The finite candidate graph goes to the existing exchange backend, and
its resulting coordinates are checked again independently using exact integer
norm, pairwise angle and rank calculations.

The `blocker_limit`, `candidate_limit`, `pool_seed` and `seconds` parameters
have the same meaning as in `search.partial_sign_orbits`. Their defaults are
respectively the source size, 20,000, no random seed, and one second. The
receipt records truncation, the tested words, and the common denominator. A
failed or timed-out search is not an upper bound or an optimality claim. The
current explicit graph and Gram work are quadratic in materialized candidates;
the existing resource limits apply, so this is not a compressed backend for
million-point constructions. It supports only rational isometries of a fixed
positive diagonal metric, not general nonlinear maps or algebraic-field
coordinates. Other construction operators remain necessary.

Run `python -m unittest discovery.tests.test_rational_images` for small exact
controls, including noncommuting words, weighted metrics and an invalid
isometry. The 14D Ganzhinov configuration was also used in a local zero-blocker
smoke test: it returned its original 1,932 points with an exact proof and no
improvement. Neither test is a new record.
