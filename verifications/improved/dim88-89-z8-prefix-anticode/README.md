# Twelve punctured Z8 fibers in dimensions 88 and 89

**12,884,901,888 kissing points in 88D**, and the same bound in 89D by appending a zero coordinate. The 88D gain is 3,294,110,112 (+34.34659%) over the preceding repository bound. The 89D gain is 194,746,184 (+1.53462%). No unrestricted optimality or exhaustive worldwide-priority claim is made.

Retain prefixes (0,b), b=0,1,2,3, in each of three crosspolytope alphabets for the length-24 rank-12 Z8 parent. Delete the first two four-coordinate blocks. All resulting points have squared norm 88. Distinct points in the same fiber have inner product at most 40; distinct same-color fibers at most 44; different colors at most 22*2=44. Thus normalization by sqrt(88) gives a kissing configuration. These strict bounds below 88 also prove injectivity and disjointness. Every prefix fiber has exactly 8^10 elements.

## Portable replay

Run `python verify.py` (Python 3.10+, NumPy). This first replays the complete inherited Z8 minimum proof, then the independently written integer alphabet, rank, fiber-count and pair-type checker. There is no network or other-package dependency. No 12.9-billion-point enumeration is required: the finite proof exhausts the algebraic pair types.

Files and scripts:
- `verify.py`: full replay driver.
- `check.py`: independent exact finite geometry and counts.
- `parent_independent.py`: unchanged inherited parent proof implementation; its unused historical sparse-layer routines are retained for provenance. Invoke via verify.py, not standalone.
- `z8_generator.json`: frozen exact parent generator.
- `parent_receipt.json`: original parent replay report.
- `claims.json`, `final.py`: current registry values and driver.
- `AUDIT.md`: complete proof, source identity and comparison limits.
- `result.json`: independent finite-check receipt, whose parent-minimum conditional status is discharged by verify.py.
- `parent_replay.log`: fresh inherited parent replay completion.

The generator and parent proof are copied unchanged from [the prior three-alphabet package](../dim93-96-three-crosspolytope-alphabets/three_basis_independent_20261117/). The full replay reuses that independently authored parent algorithm, not a newly authored second parent-minimum implementation. The source parent is attributed to Duursma, Greferath, Litsyn and Schmidt (2001), DOI [10.1109/18.923742](https://doi.org/10.1109/18.923742).

The external-source comparison in the main table remains the derived Gamma72 plus tail bound, not an independently tabulated 88D or 89D record. This package proves the explicit construction and its gain over the repository, not worldwide priority.
