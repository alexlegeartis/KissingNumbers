# 85D: fill vacant contacts on already used phases

The bound is **6,294,162,634**, adding 1,385,066 contacts to the
[previous Barnes–Wall-type tail construction](../dim82-85-bw16-tail-sections/)
without deleting any of its points. It still uses Nebe's Γ72 core.

That construction chooses 384 of 453 antipodal tails and gives each
a different horizontal phase. We now use all 69 remaining rays.
A phase can serve more than one tail if their horizontal points
are disjoint: two different points in the same coset differ by a
nonzero Γ72 vector, with squared length at least eight. If the
horizontal points coincide, however, the pair is forbidden even for
orthogonal tails. This is enforced row by row, not inferred from
phase labels alone.

Each new ray uses only 4+4+4 horizontal profiles. In a different
phase of the same colour, the only potentially short horizontal
difference is one occupied norm-eight Leech block; the other two
blocks must coincide. This cannot couple with a 4+8 or 6+6 profile,
so only already selected triples can obstruct the new rows.
Different colours are safe by depth six. The concrete extra ray
and phase order is pinned by `extra.json` and replays exactly.

The 86,558 additional triple rows contribute 16 points each
(eight horizontal signs and two tail orientations). Include also
the 138 pure contacts from the same 69 remaining antipodal rays,
compatible by the exact full tail Gram. Run
`OPENBLAS_NUM_THREADS=1 python3 verify.py`; this replays the parent
using Python, NumPy, SymPy and a little-endian C++17 compiler.
