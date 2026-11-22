# Rotated five- and three-plane layers in dimensions 86, 87 and 90

| Dimension | Previous bound | New lower bound | Gain |
|---:|---:|---:|---:|
| 86 | 8,008,545,838 | **8,008,552,282** | 6,444 |
| 87 | 8,160,981,972 | **8,160,988,416** | 6,444 |
| 90 | 21,944,100,474 | **21,944,111,386** | 10,912 |

Keep the QR48 and Pless layers and the square supports in the actual
[86–87D](../dim86-87-remarked-square/) and [90D](../dim90-91-a2-rotated-square/)
predecessors. Ordinary triangle planes use the same angular frame as
in the [88D construction](../dim88-89-signed-row-designs/).

On each of the stored five-plane supports use unit entries along the
square-frame axis 30 degrees from the nearest Pless direction, with 16
even sign patterns. Different supports share at most two planes, and
none contains an old square support, so cross dots with squares are
at most 3/√2. Pless cross dots are at most 5√3 cos 30° = 15/2, exactly
at the required threshold. QR cross dots are at most 5. The five-plane
layer removes the four signed D roots on each of the used axis pairs.

For each triple missed by the square support design and not contained
in a five-plane support, take all eight signs on the unit direction
tilted 15 degrees from the selected axis, provided all three axis pairs
are already paid for by the five-plane replacement. The omitted square
triples form a Steiner triple system, so distinct triples intersect at
most once, and any two signs on one triple give dot at most 1. Triple
dots against five-plane points are at most 2 cos 15°, and against
squares at most 2 cos 30°; no new D roots need to be removed because
the mixed-axis bound is cos 15° + sin 15° = √(3/2). The dense and sparse
dots are at most 3 and 3√3.

In 86–87D there are 543 five-plane and 90 three-plane supports,
and 2,964 retired roots. In 90D the counts are 864, 131 and 3,960.
The exact replay `python verify.py` imports and pins the actual old
support generators and validates the finite support, phase and angular
certificates. No optimizer or input download is used.
