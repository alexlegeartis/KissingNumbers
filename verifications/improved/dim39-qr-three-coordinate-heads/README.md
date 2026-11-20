# 39D: 783,468 kissing points with three-coordinate heads

This explicit construction gains **8,496** over the preceding
[774,972-point five-head code](../dim39-qr-pentapolar/). It is a lower bound,
not a proof of optimality.

Three coordinates carry a rational head and the other 36 carry signs of
magnitude 2. Shortening the length-39 binary QR kernel on the three
logical marks yields 4,096 words of minimum Hamming distance 12.
Eight head positions carry ten cosets each, and 4 further positions
carry one coset each. These 84 complete fibres give **344,064** dense
points. The locations of the physical marks need not coincide with the
logical shortening coordinates.

A body-coordinate permutation and a single fixed sign flip align some small support restrictions with dual-code relations. This lets the exact verifier retain sparse points that a magnitude-only bound would discard.

The retained sparse layer consists of 436,432 signed octads,
2,908 roots and 64 stars. The verifier removes all
conflicting signings and checks exact cross-layer angular inequalities.

Run `python verify.py` with Python 3 from this directory. Hamming distances,
dot products and squared angle tests are checked with Python integers.
No numerical optimization tolerance is used in accepting the witness.

Starting binary code: [Yang et al., arXiv:2609.35051v1](https://arxiv.org/html/2609.35051v1),
Appendix B, Tables 11–12. Original signed support data:
[Lin, arXiv:2611.03682](https://arxiv.org/abs/2611.03682), CC BY 4.0, and its
[companion repository](https://github.com/QH-Lin/kissing-numbers-33-39).
The original support bytes are gzipped; the repository's seven prior
support changes and the retained star data are listed separately.
