# Rephased sparse banks and off-axis DG donors, 66D–68D

| Dimension | Previous repository bound | New verified bound | Increase |
|---:|---:|---:|---:|
| 66 | 434,387,366 | **434,471,332** | 83,966 |
| 67 | 454,068,488 | **454,416,648** | 348,160 |
| 68 | 485,390,106 | **486,029,082** | 638,976 |

The [preceding D4 donor construction](../dim66-68-dg-d4-tail-donors/)
paid for the new dense layers by omitting 3, 170 and 312 tagged sparse
banks. Instead reflect selected coordinates separately for each bank:
a 16-bit mask switches every point's signs on those coordinates without
changing any within-bank angle. On each signed overlap between two
banks the two changes to the guard parity must agree; on each sparse
bank–atlas guard the parity must remain fixed. The witness masks solve
those simultaneous constraints over F2. They also preserve the zero-tail
banks' even parity and all tagged parity requirements of the
[older weighted triple donors](../dim65-68-dg-diagonal-triples/).

In 67–68D the donor layers are unchanged and all 170 or 312 paid banks
are restored. In 66D use four donor layers instead of two: the tails
are ±(3.1,−3.1) and ±(5.1,2), with exact cubic labels and disjoint
quadratic translates recorded in selection.json. All sparse banks
are again retained. The exact count pays for 94 atlas pieces and two
axis poles in 66D. No other choices are removed.

Run `python3 -B verify.py` without assertions disabled. It requires
a BMI2-capable x86-64 CPU and a C++20 compiler. The proof replays the
published triple geometry, checks all 31,512 shared overlap masks and
the atlas guards, then verifies every new angle and exact count.
clearance.cpp exhaustively checks the relevant binary distances,
and translates.json contains the literal sign translations.

Optional: compile and run local_ceiling.cpp with a C++20 compiler.
It computes the exact Walsh spectrum of the graph of weight-16
differences in RM(2,6): degree 2,604 and least eigenvalue −84.
The ratio bound is 131,072, attained by one full DG layer. This only
limits a single cubic and fixed tail in the stated distance regime.
