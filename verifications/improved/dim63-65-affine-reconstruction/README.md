# Reconstructed 63–65D affine sign-subcode extension

Starting with published KissingNumbers commit 69649636's global64 certificate,
preserve all layer cardinalities and add 33 smaller affine sign classes, changing 81
previously untouched parent Hamming classes without changing their sizes.
The exact gain is 38,592 points in dimensions64 and65:
335,749,696 and335,817,922. The best of the 72 tested octad constructions has184,415,700 points,
21,712 above the preceding63D bound.

`verify.py` checks every changed sign code and all relevant support-pair
interactions exactly; `section63.py` explicitly enumerates affine sign words.
Base geometry is inherited from the repository's full independent parent replay;
the64D certificate and65D owner lift must be replayed with the original package.
All source paths resolve within this repository.
The complete parent proof package is `verifications/improved/dim61-67-compatible-sign-families`.
Independent exact audit passed; see `audit/README.md` and its reports. This package does not claim global optimality.

Discovery: in `discovery/`, run `prepare.py`, compile the repository's `subcode67_collect.cpp`,
then use regenerated `codes/candidates.bin`, and run `joint.py`.
The frozen discovery pool and context are included, so `discovery/joint.py`
can replay selection immediately; its output stays in the discovery folder.
Collector finds35 candidates, total38,848 points; compatible selection yields38,592.
No claim of optimality within all possible sign constructions.

## Explicit 63D construction convention
The 63D result is a section of a relabeled variant of the 64D construction,
not necessarily the literal section of its frozen coordinate representative.
Independently relabel the quadruple-support SQS so its labels 0–7 occupy the
selected octad. `section63.json` gives the permutation and `section63.py`
enumerates its 102,900 section points. This relabeling preserves every
cross-layer inequality: quadruple vectors have amplitude 4 on four positions,
so their inner products with amplitude-2 weight-16 vectors, dense amplitude-1
vectors, and amplitude-8 axes are at most 32, 16, and 32 respectively.
Untouched full Hamming classes are independently labeled with their section
positions first, as in the original section proof. Explicitly constrained
classes retain their certified check matrices and are counted by enumeration.

## Reproduce
From the repository root, using Python and NumPy 2.x:
```
python verifications/improved/dim63-65-affine-reconstruction/verify.py
python verifications/improved/dim63-65-affine-reconstruction/section63.py
python verifications/improved/dim63-65-affine-reconstruction/audit/audit.py
python verifications/improved/dim63-65-affine-reconstruction/audit/section_audit.py
```
For the unchanged base layers and 65D owner lift, also run the original
`verifications/improved/dim61-67-compatible-sign-families/verify.py`.

`final.py` derives registry counts from fresh sign-word enumeration, not cached receipts. `discovery/utils.py` provides support and code utilities for the optional search.
