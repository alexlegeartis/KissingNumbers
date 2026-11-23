# 81D: D9 roots and an even octad on Γ72

**6,251,562,978**, adding **39,168** to the
[pinned Γ72/P9a construction](../dim81-84-antipodal-tail-codes/).
Keep every old point. New tails consist of the 144 signed (D_9)
roots, plus the 128 even signs on a fixed set of eight coordinates.
The roots split into nine round-robin matchings; the octad signs split
into eight cosets of the affine Reed–Muller code RM(1,3). Each of
these 17 sets has pairwise nonpositive dot products after normalizing
to squared norm four.

Use a distinct nonzero B-only phase for each set. Its 144 horizontal
vectors are the norm-four contacts in three Leech frames, with the
same exact separation argument as the
[82–83D norm-four layer](../dim82-83-gamma72-norm4-tail-weave/).
The two exact gates check all new–new pairs and all old tail
directions: B-only phases separate by squared distance four, and
different mixed colors by squared distance six. The Γ72 equator
clears by the norm-four phase depth, and old pure tails are twice the
old mixed tails. No discovery solver is needed.

Run `OPENBLAS_NUM_THREADS=1 python verify.py` from this directory
with Python 3, NumPy, SymPy and a C++17 compiler. The default also
replays the pinned parent and asserts the literal frames; use
`--incremental` to check the new tail geometry separately.

Within the norm-eight signed-octad tail chart the 128 octads are
optimal, even allowing their missing coordinates to vary. A signed
octad is an edge of the nine-cube along its missing coordinate.
Distance compatibility says that the chosen edges form an induced
matching: a different omitted coordinate leaves seven shared bits
and requires two disagreements, as do the eight shared bits of two
parallel edges. If (S) is the set of matching
endpoints, each has eight neighbors outside (S). A point outside
(S) cannot have all nine neighbors in (S): the partner of any
one of those neighbors would then touch a second one. It has at most
eight. Counting boundary edges gives (|S|le 256), hence at most
128 octads. This is not an upper bound for mixed tail types or for
the 81D kissing number.
