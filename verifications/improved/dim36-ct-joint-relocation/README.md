# Joint Coxeter–Todd trine relocation in 36 dimensions

Begin at the [502,929-point two-trine construction](../dim36-ct-double-trine/). Select 2,142 additional norm-twelve root triples in three Coxeter–Todd blocks, pricing the exact union of the blockers rather than charging each triple separately. Some groups of new triples share one old heavy vector but block its small-root trine in the original block. Relocate a complete small trine to the other, previously empty block in 203 such slots. Subsequent exact checks add another 31 points in 22 vacated slots.

The result inserts 2,142 root triples, deletes 1,648 old two-block points and 48 old root triples, and restores 640 two-block points. Thus the lower bound is **504,015**, up 1,086 from the common ancestor and 467 from the [previous pocket-trade record](../dim36-ct-pocket-trade/). All added points have squared norm twelve.

The integer replay pins the common verified ancestor, reconstructs the original population and explicitly checks each new Gram product. The literal witness identifies removals, relocations and remaining partial phases. There is no numerical tolerance or optimization in the verifier. Run `OPENBLAS_NUM_THREADS=1 python3 verify.py` with NumPy. This is a lower bound, not an optimality claim.

For a separate literal full-coordinate Gram replay, run `OPENBLAS_NUM_THREADS=1 python3 audit.py`.
