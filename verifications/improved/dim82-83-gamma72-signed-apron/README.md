# Larger signed norm-4+8 packets in 82D and 83D

| Dimension | Previous repository bound | New verified bound | Increase |
|---:|---:|---:|---:|
| 82 | 6,269,084,390 | **6,269,087,562** | 3,172 |
| 83 | 6,278,756,940 | **6,278,767,980** | 11,040 |

This refines the [independent signed tails](../dim82-83-gamma72-signed-packets/)
of the Γ72 construction, keeping the tail code, both phase labels and the
4+4+4 and 6+6 horizontal profiles. The 4+8 choices are now allowed to keep
the sign of the norm-eight block independently, while preserving both
signs in the norm-four block. A selected mask bit thus contributes two,
not four, points.

The only unfavorable different-phase case has a positive tail dot and
a single raw norm-eight difference in the horizontal part. A single
changed block cannot connect different norm profiles in these packets:
to preserve the norm of the whole horizontal vector, the changed blocks
must have equal squared norms. When the norm-four block changes, identical
signed norm-eight vectors produce the only overlap to be checked.
When the norm-eight block changes, the risky unsplit local inner product
is +4; compare equal signs at +4 and opposite signs at -4 on canonical
representatives. With equal phases check literal coincidence. Antipodal
tail signs have squared distance eight and require no horizontal exclusion.

Run `OPENBLAS_NUM_THREADS=1 python3 -B verify.py` with NumPy and SymPy.
The script replays the old packet proof first, then checks the new bit
masks and exact resulting counts. It expects the sibling verification
directories; pass their parent explicitly as the first argument if
running elsewhere. The construction depends on the same named Γ72 and
small-dimensional tail results as the preceding package.
