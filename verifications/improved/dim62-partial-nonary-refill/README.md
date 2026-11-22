# Partial nonary signed refill in 62 dimensions

**126,571,579 kissing points**, an increase of **943** over the
[preceding weighted-distance proof](../dim62-coordinate-weighted-distance/).
We keep all dense points, all weighted omitted-fibre spikes and all
other sparse classes. The certificate inserts 1,914 nine-support
vectors on 60 new supports and retires 971 identified signs of the
old nine-support layer. No other deletions are needed.

Work at squared norm 64 and dot-product cap 32. The new points have
amplitude 8/3 on nine coordinates. Sign sets in the same cell are
subsets of a length-nine, distance-three sign code; the literal
minimum is checked rather than inherited from its discovery. Signed
overlap at most four is sufficient for distinct cells and old
nonary rays, since 4(8/3)²<32. For every overlap above four the
verifier checks every relevant sign pair, not just supports. It
checks the full old-deletion union, with duplicates removed.

The weight-sixteen classes have amplitude two; intersection ≤6
gives dot at most 32. The amplitude-four classes are affine
quadruples with zero coordinate XOR. Distinct pairwise XORs in
each new support forbid containing any such block, so their
dot is also at most 32. Retained poles have amplitude eight,
giving dot at most 64/3. The ordinary dense latitudes have
body amplitude one and two tail coordinates of magnitude
at most two, so their dot is at most (8/3)·11. The omitted-fibre
spikes have ordinary absolute values b or (9/2)b and zero
tails, where b²=256/317. Their new-point dot is at most
(8/3)·(25/2)b<32. These are universal envelopes, independent of
the fixed phases or which duplicate spike words were kept.

Run `OPENBLAS_NUM_THREADS=1 python verify.py` in this directory.
It pins and replays the parent and rebuilds the complete lists
of its old sixteen- and nine-support classes from the linked
proofs. The selected signs and paid individual blockers are
literal in [winner.json](winner.json). The search optimizer
is not needed for replay. This is a construction lower
bound, not an optimum claim.
