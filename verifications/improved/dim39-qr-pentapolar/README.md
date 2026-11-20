# 39D: 774,876 kissing points with five head positions

This is an explicit lower bound, **5,024** above the preceding repository
witness of 769,852 and **5,960** above the named 768,916-point Lin construction.
It is not an optimality claim.

Mark two coordinate axes, leaving a 37-coordinate sign body. Shorten the
binary QR kernel from 39 to 37 coordinates to obtain 8,192 words, with
minimum Hamming distance 12. In the two marked axes, use five rational heads.
The first four heads each get ten cosets of the shortened kernel. The fifth
head, farther from the origin, gets one additional coset; its relative
distance to the four groups is checked exactly. Each dense vector has
body coordinates ±2 and the chosen rational two-coordinate head.

This gives 335,872 dense points. Keep the previous 3,424 octad supports,
but remove the even signings with both marked entries positive. Exactly
70 octads contain both marked axes and each loses 32 points, leaving
436,032 octad points. The 2,908 roots and 64 stars complete the witness.

Run `python verify.py` with Python 3. The replay checks all within-fibre
and cross-fibre dense angles via exact Hamming minima, checks the sparse
support and root/star inequalities, and checks every sparse-to-head gate.
All angular tests use signed integer products and squared norms;
there is no tolerance or floating-point acceptance.

## Sources

- [Yang et al., Large-Scale Autonomous Discovery of Kissing Number Constructions](https://arxiv.org/html/2609.35051v1), Appendix B, Tables 11–12, for the starting QR generator matrix and representatives.
- [Qihong Lin, New lower bounds for kissing numbers in dimensions 32–35 and 37–39](https://arxiv.org/abs/2611.03682), CC BY 4.0, and the [accompanying data repository](https://github.com/QH-Lin/kissing-numbers-33-39) for the 3,424 original sparse supports. The compressed file contains Lin's original unaltered support list, with the retained repository edits given separately in `base_certificate.json`.
