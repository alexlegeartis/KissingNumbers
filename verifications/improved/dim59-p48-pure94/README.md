# 59D: sparse antipodal pure-pole replacement

**56,969,157 = 56,969,091 + 66.** Keep all old equator and mixed
points. Replace the 28 pure poles surviving the preceding weighted
cube interface with 94 literal poles on 47 antipodal lines.

For each nonzero integer row y in the witness add the two unit
tail directions ±y/||y||, with zero 48-dimensional body. Every
pure-pure comparison is checked by 4⟨y,z⟩²≤||y||²||z||². Their
interface to the long 48-dimensional owner stock needs no scanning:
every old lift puts exactly amplitude 1/2 in the tail, so the
Cauchy bound ≤1/2 applies to every unit pure pole. The equator
has zero tail. For the inherited tetrahedral signed triples,
the predecessor's 3-support, sqrt(5/24) coefficient and each
literal pole's support L1 norm certify all pairs. The inherited
unequal-coordinate cube uses weights (50×8, 49×3); the
corresponding weighted L1 bound certifies all its signs.
The 32 smaller mixed tails are checked individually on their
pinned rational coordinates. The verifier reconstructs exactly
the 28 previous survivors; no claimed saving on pure-pole
deletion not actually checked is used.

Run with NumPy:
  OPENBLAS_NUM_THREADS=1 python verify.py --parent WEIGHTED_CUBE_DIR --owners POINTWISE_P48_DIR --body GLOBAL_OWNER_REFILL_DIR
