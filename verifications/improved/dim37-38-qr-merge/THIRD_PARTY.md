# Third-party inputs
- **Qihong Lin,** New lower bounds for kissing numbers in dimensions
  32–35 and 37–39, [arXiv:2611.03682v1](https://arxiv.org/html/2611.03682v1).
  The unchanged input `supports_3088.txt` is his published code
  `certificates/d38/supports_3088.txt` at
  [pinned revision 0c4e2bd](https://github.com/QH-Lin/kissing-numbers-33-39/blob/0c4e2bd486061b3fd55e128ac3497bfb12b73f07/certificates/d38/supports_3088.txt).
  It carries [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/);
  see the [repository license](https://github.com/QH-Lin/kissing-numbers-33-39/blob/0c4e2bd486061b3fd55e128ac3497bfb12b73f07/LICENSE).
  The SHA-256 matches that published in Lin's Table 4.
- **S. Yang, R. Zhao, J. Wu et al.,** Large-Scale Autonomous Discovery of
  Kissing Number Constructions,
  [arXiv:2609.35051v1](https://arxiv.org/html/2609.35051v1), Appendix B,
  Tables 11–12, [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
  The length-39 generator and coset representatives are transcribed
  unchanged into JSON. The new coordinate compression is not part of
  that source.
  The length-38 parent masks in `dense_parent38.json` are from the
  same two tables, with the same license and attribution.
- The 36-coordinate numerical support used in 37D is the construction
  tabulated by A. E. Brouwer as
  [A(36,8,8) ≥ 2742](https://aeb.win.tue.nl/codes/Andw.html).
  Instead of redistributing the original [list](https://aeb.win.tue.nl/codes/cwc/d8/a36.8.8.2742H),
  `supports37.py` gives an explicit orbit formula and confirms the
  sorted numerical set equals that list by digest.
