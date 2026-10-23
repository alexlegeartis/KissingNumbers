# 43D certification status

## Established: 2,878,155 points

`certify43.py` certifies the recovered `data/old43_seed.npz`. A fresh run on
23 October 2026 passed:

- Exact reconstruction of the pinned P48n embedding, its full rank-40
  complement and all discriminant classes and gluing witnesses.
- Fresh exhaustive integer enumeration of every one of the seed's thirteen
  large body shells: twelve at squared radius 16/3 and one at radius six.
- Exact minimal-shell 11-design identities establishing 9,180 bodies in each
  of 41 families and 495 bodies in each of 59 families.
- Exact head rank three and every one of the 113*114/2 = 6,441 internal and
  cross-family angle conditions, including coincident bodies in reused shells.
  Integer products use Python arbitrary-size integers; the head denominator
  is 10^10, which would overflow int64.

The resulting accounting is

    1,084,590 + 1,387,980 + 41*9,180 + 59*495 = 2,878,155.

The enclosure used by the integer enumerator is regenerated and justified
from an exact bound on its LDL rounding error. Floating-point Cholesky only
proposes an enclosure; exact arithmetic proves that it contains the entire
target shell. No optimization solver or inherited modular population formula
is needed for this proof. Published assumptions remain the extremal P48n
minimum six, its minimal-shell size 52,416,000 and the 11-design property.

The recorded run is `certification/recovered43.json`; regenerate it using
`python certify43.py --jobs 4 --output NEW_DIRECTORY`. The default consistency
verifier still labels its own population check as inherited. Recorded fresh
counts are evidence of the separate full run, not a new enumeration on each
default invocation.

## Unresolved: 2,911,990 points

The supplied writeup describes a 245-family baseline:

    1,084,590 + 1,387,380 + 45*9,180 + 50*495 + 2,163 = 2,911,983.

It translates three body classes together on a triangular face, repacks the
auxiliary heads, then adds seven points with nonzero bodies using a finer
32+11 decomposition. These are not seven pure-head additions to the recovered
seed. Their exact maxima against all original families must be established.

The writeup names these external inputs, whose coordinates it does not print:

- `kissing41-46_exact_certificates.zip`, containing the 245-family baseline.
- `kissing43_extension_2911990.zip`, containing
  `packed_extension_2911990.npz` and `verify_extension.py`.

Neither named archive is supplied in the two research ZIPs. The head selection
for the 45 large auxiliary families, all smaller baseline layers, and the
seven new points cannot be inferred from the population totals. Rounded
maximum inner products in a writeup do not supply an exact certificate.

Until those witnesses are restored or a new construction is independently
generated, the 2,911,990 figure remains a documented working claim and is
excluded from the results registry. The older seed is not substituted for it.
