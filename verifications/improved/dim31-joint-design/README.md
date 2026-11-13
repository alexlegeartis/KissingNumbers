# 31D: joint design and exterior deletion optimization

**251,396 points**, one more than the predecessor 251,395 and 142 more than PackingStar commit 2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6. This is a dated pinned comparison, not an exhaustive global-priority claim.

Instead of rotating one fixed Leech–Sloane module, choose its 70 transversal quadruples jointly with the surrounding deletion cost. The 280 triple equations select a Steiner quadruple system; 560 possible blocks each contribute 16 signed points. The fixed part has 1,444 points. At the preceding best orientation, the full new module requires 140 exterior deletions rather than 142. Its union with the predecessor yields 1,004 additions and 1,003 deletions.

`python verify.py` reconstructs the 2,564-point module from the exact-cover design and two affine binary glue codes, checks the rational isometry, independently checks every added point against every final point with int64 arithmetic, and checks all 31,599,848,710 unordered final pairs using bounded exact integer products in binary64. All source coordinates and exchange data are included. No optimizer result is assumed by the proof.

Search source and numerical MILP logs are in the research checkpoint. Solver optimality is not used as an exact family-upper-bound claim. The earlier proof packages remain historical witnesses.

## Files and reproduction

- `final.py` registers the count read from the full coordinate archive.
- `module_recipe.py` reconstructs the exact Leech–Sloane module from the selected design and binary glue codes.
- `verify_geometry.py` performs complete bounded-integer all-pairs verification, grouping antipodal pairs only when both signs occur.
- `verify.py` verifies the module orientation and exchange and invokes the full geometry check. Run from any working directory with Python 3 and NumPy installed.
- `data/` contains both full source and final coordinate archives, the exact design and orientation, the exchange, and verification reports. Reports are receipts, not proof assumptions.

The coordinate archives store twice the displayed norm-32 coordinates: squared norm 128 and inner-product limit 64, with denominator 2. All proof arithmetic is integer or exactly represented bounded integer arithmetic; no angle tolerance is used.
