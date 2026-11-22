# Dimension 37: 544,538-point QR code with a 12-coordinate phase refill

Run `python verify.py`. The standard-library-only replay verifies the 544,474-point
predecessor and all new interactions exactly. This raises the predecessor by 64.

Remove just the 384 parity phases on the three 8-supports of the ordinary
12-flat U=12..23 (zero-indexed). Replace them by 448 signed weight-eight
vectors, supported on eight different masks. All nonzero coordinates are ±2.
The 14 phase orbits in `win.txt` give a compact reconstruction: the first
line is a basis of a six-dimensional binary sign-flip group on U; subsequent
lines specify support, initial negative mask, and orbit size. The original
QR and radical layers are retained outside this replacement.

Every remaining exterior ordinary octad and radical six-support meets U
in at most four positions. Thus the raw cross product against the new
weight-eight vectors is at most 16. Root products are at most 16; poles
are orthogonal; dense products are at most 32, with dense squared norm
at least 145. New inner products use only integer signed overlaps,
exhaustively checked. The count is 544,474 − 384 + 448 = 544,538.

`search.cpp` is the optional seeded exploratory search (g++ -O3, C++17).
Its rank-six phase groups are tested by coset distances on support intersections.
The search is not used by the proof. The two predecessor verifiers are copied
from the existing QR and radical packages (the radical replay merely has local
import/data paths adapted).
