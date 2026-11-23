# Coxeter–Todd whole-cell reallocation in 36 dimensions

The three-block construction has a layer with squared block norms `(4,4,4)` and another with `(8,4,0)` up to permutation. First hold the `(0,6,6)` and axial layers fixed. In each placement of an `(8,4,0)` vector, the short-coset rules identify a small-root trine matched to the heavy vector. This makes it possible to price a whole replacement of `(4,4,4)` label cells by recomputing which trine phases each heavy vector can carry, rather than pricing isolated root points against the original heavy layer.

Starting from the [502,929-point ancestor](../dim36-ct-double-trine/), the literal witness adds 3,077 root triples, removes 210 old triples and 2,396 old two-block points, and places 965 reassigned two-block points. The result is **504,365**, a gain of **202** over the [504,037-point min-cut construction](../dim36-ct-joint-relocation/).

The replays reconstruct the retained and added points and verify separation with exact integer Gram products. The declared removals may include additional points beyond the currently visible blockers, and the counts debit all of them. No optimization or assumed optimality is used in either replay. Run `OPENBLAS_NUM_THREADS=1 python3 verify.py` or, for the separate literal-coordinate audit, `OPENBLAS_NUM_THREADS=1 python3 audit.py`. Requires NumPy.
