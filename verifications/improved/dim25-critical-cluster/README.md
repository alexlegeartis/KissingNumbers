# Critical-cluster refinement: 197,586 points in dimension 25

9 November 2026. This exact certificate improves the preceding repository
construction from **197,585 to 197,586** points, and is two above 197,584.
It consists of:

- **194,561** retained vectors of the 196,560-point Leech shell;
- **2,778** points on 926 independently rotated circles;
- **247** rational general points.

Every point has squared norm **32**, and every distinct-point inner product
is at most **16**. Thus normalization gives a kissing configuration in R^25.
The preceding [197,585-point](../dim25-point-split/) and
[197,584-point](../dim25-collective-deformation/) packages are preserved.
This is a constructive lower bound, not an optimality result or a fresh
worldwide-priority claim. The manuscript is unchanged.

## Reproduce the complete proof

Python 3.11 or newer and NumPy suffice. The tested environment is Python
3.12.14 with NumPy 2.3.5. From this directory:

```sh
python -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python final.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verify.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python self_test.py
```

Pin compute jobs to an allowed CPU while leaving at least one logical CPU
free, as required by the repository's AGENTS.md. A complete replay took
about ten seconds on one CPU on the publication host. No optimizer,
network access, stored shell, construction module, or floating tolerance is
needed. The package is self-contained and can be copied out of the repository.

`verify.py` checks the complete SHA-256 manifest and copies only the checker
and certificate to a temporary directory before replaying the proof. It
does not rewrite evidence. `--output /tmp/fresh25.json` optionally saves a
fresh receipt outside the package. The checker can also be run directly:

```sh
OPENBLAS_NUM_THREADS=1 python source/verify_general.py \
  source/certificate25.json --output /tmp/fresh25.json
```

`data/verification.json` and `data/verify.log` record a fresh publication
replay. `self_test.py` tests exact-sign arithmetic and rejects seven
corrupted witnesses, including a same-norm duplicate and a retained shell
blocker. `final.py` derives the count from the actual circle and deletion
lists for the repository registry; it is not the geometric proof.

## Proof scope

The standalone checker rebuilds the binary Golay code from its polynomial,
checks its dimension, weight distribution and self-orthogonality, and
enumerates all three Leech-shell families. [PROOF.md](PROOF.md) includes the
complete elementary shell/shell argument, rather than assuming a lattice
minimum. Circle/shell pairs follow from checked endpoint membership and
convexity. Different same-class circles use a checked centre bound.

The remaining comparisons use outward integer intervals with explicitly
bounded integer matrix products, followed by exact signs in
Q(sqrt(2),sqrt(3)) whenever the intervals are inconclusive. All **4,573,800**
exceptional pairs are accounted for, including **1,658,754** same-class
circle pairs and **2,915,046** interval-screened pairs. There are **3,292**
exact exceptional fallbacks. All **48,056,567** rational-point/retained-shell
pairs are checked, with **12,492** exact fallbacks. Norms, count and pair
bounds together also prove uniqueness.

## Construction and exactification

A critical cluster is freed from the circle parameterization and replaced
by general points, with additional shell deletions. Relative to the
preceding published count, this certificate has 45 fewer circle triples,
144 more general points and eight fewer retained shell points, a net gain
of one. These are category counts, not a claim that all remaining points
have unchanged coordinates.

Inherited contacts must be preserved exactly. Ordinary rational rounding
of an intermediate candidate failed. Eighteen freed inherited-circle
points were instead represented by rational points on their exact conics;
general point 231 was projected back to original circle point 2609
(displacement approximately 6.923e-8 in the numerical search coordinates)
and rationally parameterized on that conic as well. The resulting 19 conic
points and all other general points have exact rational norm 32. Acceptance
comes from replaying every pair of the final witness, not from the numerical
displacement, search residuals, or exactification procedure. No rejected
candidate is included in this package.

`source/certificate25.json` is the exact published witness. Its SHA-256 is
`5270afec106a2336179fd8c1d21fd5b9e30da100a8d4f1e6ec4ca3ec8b161299`.
The standalone checker is byte-for-byte unchanged from the preceding
package (SHA-256
`4becb67a85b17ae419570d9d10da4c931ccc5dbd758573756f48d537850d4a86`).
Separate full replays passed before publication. This is independent
reconstruction and replay, not a claim of two independently implemented
checkers.

Each class lists a 24-dimensional shift and shell endpoints (`owners`).
`local_owners` and `new_owners` list additional deletions. Three circle
points are reconstructed per owner, using the exact phases specified in
the checker and rational half-angle numerators with their common positive
denominator. Each `small` entry supplies 25 integer numerators and a
positive integer denominator. Large integers are decimal strings.

The witness retains its historical candidate description and exactification
metadata. These fields and its source, parent and numerical-seed hashes are
provenance only, not replay dependencies. The checker's `record_claim: false`
means it proves geometry rather than assessing literature priority. Its
comparison field remains the historical 197,584 baseline; the gain over
published 197,585 is one, as recorded in `data/provenance.json`.
