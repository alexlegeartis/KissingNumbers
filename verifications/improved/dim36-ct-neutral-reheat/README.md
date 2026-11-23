# 36D: coordinated cell perturbations and exact recovery

The compact Coxeter–Todd witness gives **506,628** points of squared
norm 12, with distinct pair inner products at most 6. It improves the
preceding 506,376-point construction by **252**. From the pinned
502,929-point stock it adds 14,597 triple-root points, removes 7,496
old triple-root and 5,850 old two-block points, and inserts 2,448
reassigned two-block points. Temporary losses during the search opened
up combinations of whole-cell moves that subsequently recouped
more than those losses; the proof needs only the final literal witness.

```sh
OPENBLAS_NUM_THREADS=1 python3 verify.py
```

For the independent exact 36-coordinate delta replay, first replay the
independent stock in [the preceding package](../dim36-ct-multicell-fibres/).
Then run:

```sh
OPENBLAS_NUM_THREADS=1 python3 verify_independent_delta.py ../dim36-ct-multicell-fibres/witness.json witness.json
```

Witness SHA-256:
`5327877dfd99b826787391c9ecdb76db9f58747ae693f8540ef03f0b3b678c93`.
