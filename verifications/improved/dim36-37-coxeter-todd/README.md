# Coxeter–Todd phase codes in dimensions 36 and 37

Checkpoint: **24 October 2026**.

| Dimension | Cohn table, fetched 24 October | This construction |
|---|---:|---:|
| 36 | 484,760 | **501,633** |
| 37 | 500,008 | **505,359** |

The [dated table comparison](../../../common/cohn_20261024.json) records the
current numerical bounds. These are constructive lower bounds; no optimality
is asserted.

## Replay

From this directory, using Python 3.9+ and its standard library:

```bash
python final.py
python verify.py
```

From the repository root:

```bash
python run_all.py --repository-only --only coxeter-todd
```

`final.py` counts the fixed base, three-point phase orbits and latitude cells.
`verify.py` checks the complete SHA-256 manifest and runs `phase/verify.py` in
a temporary copy, preserving the supplied reports. The exact replay derives
the populations and checks every required angle; it uses no optimizer.

## Construction and proof scope

Take three copies of a rank-12 Coxeter–Todd lattice K with Gram matrix G.
An integer matrix M satisfies `M G M^T = 3 G` and `|det M| = 729`. Rational
exhaustive enumeration proves the short-shell populations 756, 4,032 and
20,412 at squared norms four, six and eight. It proves the minimum of K
and the minimum 12 of its sublattice L, and determines all minima and
multiplicities of the 729 quotient cosets.

The fixed part consists of **403,704** norm-12 vectors in 1,429 outer cosets.
The isometry `Phi = I - M` has order three. The certificate selects **32,643**
three-point phase orbits over 3,979 product cells, giving

`403704 + 3*32643 = 501633`.

Coset distances handle pairs whose three labels differ; exact integer Gram
products handle the remaining phase pairs. The verifier checks **1,020,306**
base pairs, **5,685,991** phase/base coset pairs and **6,563,399** phase-orbit
pairs. Distinct vectors in the same sublattice coset have squared distance
at least 12. Every point has squared norm 12, so normalization gives
inner products at most one half.

For 37D, 207 norm-eight body cosets each supply nine vectors. Appending
heights +2 and -2 gives two latitude layers, each of size 1,863:

`501633 + 2*207*9 = 505359`.

The replay checks **21,321** latitude-body coset pairs and **1,119,456**
latitude/equator coset pairs. Norms and distances prove compatibility
within and between the layers and with the equator. The
[source proof](phase/proof.txt) gives the complete argument. No parent
minimum, shell population or external code table is assumed.

## Artifacts and history

The files in `phase/` are byte-preserving copies of
`research/dim33_39_20261024/result/`, as recorded in `UPDATE_PROVENANCE.json`:

- `certificate.json`: Gram data, fixed cosets, phase orbits and latitude cosets.
- `exact_lattice.py`: exhaustive enumeration using rational branch bounds.
- `verify.py`, `verification.json`: original lower-bound verifier and recorded report.
- `proof.txt`, `README.txt`: mathematical proof and source instructions.
- `vectors36.npz`, `vectors37.npz`, `coordinates_check.json`: optional coefficient
  exports and their norm/distinctness checks.
- `verify_sections.py`, `section_obstruction.json`: optional NumPy check of
  moment identities excluding improvements from sections in dimensions 33–35.
  This additional negative result is not required for either lower bound.

Rows in the optional exports are coefficients in the supplied Gram basis.
Cartesian unit coordinates are `coefficients @ cholesky(gram) / sqrt(12)`.
The default replay does not require NumPy or these exports.

The earlier full-cell construction remains in `source/`: **501,147** in 36D,
and **501,149** after adjoining two poles in 37D. Its `verify.py` independently
enumerates the short shells and all 12,688,203 outer-coset pairs;
`certificate.json`, `verification.json`, `proof.txt`, `discovery.json`,
`README.txt`, `vectors.npz` and `coordinates_check.json` preserve its data
and reports. These files are retained history, superseded by the phase code.

`data/SHA256SUMS.txt` covers every installed file except itself.
