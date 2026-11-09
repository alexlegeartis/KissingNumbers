# Point splitting: 197,585 points in dimension 25

9 November 2026. This exact certificate improves the preceding repository
construction from **197,584 to 197,585** points. It consists of:

- **194,569** retained vectors of the 196,560-point Leech shell;
- **2,913** points on 971 independently rotated circles;
- **103** rational general points.

Every point has squared norm **32**, and every distinct-point inner product
is at most **16**. Thus normalization gives a kissing configuration in R^25.
The preceding [197,584-point construction](../dim25-collective-deformation/)
is preserved. This is a constructive lower bound, not an optimality result
or a fresh worldwide-priority claim. The manuscript is unchanged.

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
about five seconds on the publication host using one CPU. No optimizer,
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

`data/verification.json` and `data/verify.log` record the publication replay.
`self_test.py` tests exact-sign arithmetic and rejects seven corrupted
witnesses, including a same-norm duplicate and a retained shell blocker.
`final.py` derives the count from the actual circle and deletion lists for
the repository registry; it is not the geometric proof.

## Proof scope

The standalone checker rebuilds the binary Golay code from its polynomial,
checks its dimension, weight distribution and self-orthogonality, and
enumerates all three Leech-shell families. [PROOF.md](PROOF.md) includes the
complete elementary shell/shell argument, rather than assuming a lattice
minimum. Circle/shell pairs follow from checked endpoint membership and
convexity. Different same-class circles use a checked centre bound.

The remaining comparisons use outward integer intervals with explicitly
bounded integer matrix products, followed by exact signs in
Q(sqrt(2),sqrt(3)) whenever the intervals are inconclusive. All **4,546,620**
exceptional pairs are accounted for, including **1,814,481** same-class
circle pairs and **2,732,139** interval-screened pairs. There are **31** exact
exceptional fallbacks. All **20,040,607** rational-point/retained-shell pairs
are checked. Norms, count and pair bounds together also prove uniqueness.

## Construction and certificate format

An alternate 49-owner arrangement of the preceding family permits replacing
one general point by two, followed by joint relaxation of all 2,913 circle
phases and 103 general points. The number of deleted shell points stays
1,991. Rational half-angle rotations and rational stereographic projection
turn the resulting configuration into exact coordinates. Search diagnostics
are not used to accept this witness.

`source/certificate25.json` is the exact published witness. Its SHA-256 is
`c5ddce2af91688875d54f68ff570a2814d63dcdc7f5884ce49526b630d409d55`.
Two separately snapshotted and converted copies of the numerical seed gave
these same bytes and both passed the standalone checker. This is a separate
replay, not a claim of two independently implemented checkers.

Each class lists a 24-dimensional shift and shell endpoints (`owners`).
`local_owners` and `new_owners` list additional deletions. Three circle
points are reconstructed per owner, using the exact phases specified in
the checker and rational half-angle numerators with their common positive
denominator. Each `small` entry supplies 25 integer numerators and a
positive integer denominator. Large integers are decimal strings.

The certificate's historical description calls it a candidate; acceptance
comes from the independent proof, not that label. Its source, parent and
numerical-seed hashes are provenance only and are not replay dependencies.
The checker's `record_claim: false` likewise means the checker proves
geometry rather than assessing literature priority.
