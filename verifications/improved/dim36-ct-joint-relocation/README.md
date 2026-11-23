# Joint Coxeter–Todd trine relocation in 36 dimensions

Begin at the [502,929-point two-trine construction](../dim36-ct-double-trine/). Select 2,190 additional norm-twelve root triples in three Coxeter–Todd blocks, pricing the exact union of the blockers rather than charging each triple separately. Some groups of new triples share one old heavy vector but block its small-root trine in the original block. Relocate a complete small trine to the other, previously empty block in 203 such slots. Subsequent exact checks add another 31 points in 22 vacated slots.

A four-cell refinement treats each compatible new triple as a unit benefit and each blocker in the retained configuration as a shared unit cost. A source–sink min-cut selects the triples and pays each shared blocker only once. Compared with the [504,015-point intermediate](https://github.com/alexlegeartis/KissingNumbers/tree/220a01c660ebe09c1a80a4a4856ee7ced7ce0803/verifications/improved/dim36-ct-joint-relocation), it adds 48 triples and removes 26 points.

The final construction inserts 2,190 root triples, deletes 1,671 old two-block points and 51 old root triples, and restores 640 two-block points. Thus the lower bound is **504,037**, up 1,108 from the common ancestor and 489 from the [previous pocket-trade record](../dim36-ct-pocket-trade/). All added points have squared norm twelve. The cut is a scoped selection rule, not a global upper bound.

The integer replay pins the common verified ancestor, reconstructs the original population and explicitly checks each new Gram product. The literal witness identifies removals, relocations and remaining partial phases. There is no numerical tolerance or optimization in the verifier. Run `OPENBLAS_NUM_THREADS=1 python3 verify.py` with NumPy. This is a lower bound, not an optimality claim.

For a separate literal full-coordinate Gram replay, run `OPENBLAS_NUM_THREADS=1 python3 audit.py`.
