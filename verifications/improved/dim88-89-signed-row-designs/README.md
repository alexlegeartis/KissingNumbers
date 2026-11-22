# Exact square and triality layers in dimensions 88 and 89

| Dimension | Immediate predecessor | Net gain | New lower bound |
|---:|---:|---:|---:|
| 88 | 13,357,769,652 | 24,268 | 13,357,793,920 |
| 89 | 15,177,390,082 | 2,069,008 | 15,179,459,090 |

Run `python3 verify.py` (standard library, without optimization). The previous independently certified dense and sparse layers are inherited premises.

## 88D: better quadruple supports for the same geometry

In the [published 88D QR48–Pless code](https://github.com/alexlegeartis/KissingNumbers/tree/main/verifications/improved/dim88-axial-pless) replace only the 2,604 four-plane support sets of the rotated-square vectors. Do not change the plane frame, existing 128 even sign patterns per support, QR48 vectors, Pless vectors or D88 roots. From the affine Steiner quadruple system of order 8, take its Hanani expansion to order 22. Double it to order 44: above every old block take the eight even binary lifts, and above every old pair take its four-point inverse image. The checker directly verifies all triple incidences for the resulting Steiner quadruple systems. Deleting the two last points leaves 2,730 blocks on the 42 ordinary triangle planes. Shared support size is still at most two, so the square-square and all mixed checks in the predecessor apply without change. The exact appended-minus-deleted square count is `(2730−2604)×128`.

### Five-plane and three-plane layers

The square-frame axes are 60° and 150° relative to a Pless root axis.
Use the 150° axis, at least 30° from all Pless rays. Place unit entries
on each of 689 five-plane supports, with all 16 even sign phases. Distinct
supports intersect in at most two planes and none contains a square
support. Hence the mixed square dot is at most 3/√2, while the Pless dot
is at most 5·3/2, right at the threshold. Dense QR dot is at most 5.
Remove exactly 3,444 old D88 roots: the four signs on each of 861 pairs
of chosen axes. This adds 11,024 and removes 3,444, a net gain of 7,580.

The missing square triples are lifts of the derived Steiner triple system
at the deleted pair. On each of its triples, use the quaternary single-parity
code (symbol = lift and sign), retaining only triples not contained in
the five-plane stock. Use unit directions tilted 15° from the chosen axes.
The signed pair resource is used at most once, so mutual products are at
most 1. The square dot is at most 2 cos 30° and the five-plane dot at most
2 cos 15°, both within their normalized bounds. Against D88 roots the
mixed-axis sum is at most cos 15° + sin 15° = √(3/2), and roots using two
chosen axes have already been removed. The QR and Pless dot caps are
3 and 3√3. These 560 points require no extra removals.

## 89D: 24-cell row alphabet in four occupied rows

Keep the entire [published 89D 14-fiber plus variable-puncture Golay construction](https://github.com/alexlegeartis/KissingNumbers/tree/main/verifications/improved/dim88-89-golay-variable-punctures). Take 385 four-point supports from the same Steiner quadruple system of order 22, on those 22 ordinary 4D rows. Set the 89th coordinate to zero.

Per support choose one row frame: the signed axial directions (±2 on a single coordinate), the even-flat frame (all 16 signs restricted to even parity) or the odd-flat frame. Each consists of four orthogonal axes with signs, all row squared norms 4. Label its axes by the four binary pairs. Choose 64 length-four axis words with zero XOR and all sixteen sign words. Two different axis words of the same frame agree in at most two coordinates. For equal axis words, different sign words reverse at least one position, giving raw dot at most 8. Different frames have row dot product at most 2. Hence all same-support dots are at most 8; other supports meet in at most two rows and also have dot at most 8. The squared norm is 16.

The dense rows use these same three frames. The old Golay rows have two magnitudes 2. In either case each new active row contributes at most 4, so the raw mixed dot is at most 16, safe against old squared radii ≥80. The pure poles are orthogonal. Finally append the 24 ordinary D4 roots independently in each of the 22 rows, with zero head; their nonzero coordinates are two signed units. Their dot product with a new triality row is at most 2 and against an old Golay row at most 4, both within the relevant angular thresholds. Dense vectors, old Golay vectors, new four-row vectors, single-row roots and pure poles are distinct by ordinary support. Thus the exact gain is `385×3072+22×24`.


### Five-row tagged D4 roots

Add 865 patterns on five rows, each row tagged by one of the three
coordinate pairings. A tag gives four orthogonal root axes with entries
of magnitude √2 and row norm squared 4. Use the 64-word GF(4) projective
quadratic evaluation code, with distance three, and 16 even sign phases.
This gives 1,024 points per pattern and squared norm 20.

The overlap size plus number of matching tags is at most five between
new patterns and at most seven against every old weight-ten/eleven
Golay pattern. This gives the respective dot bounds 10 and 14√2.
Each new support avoids containing a four-row triality support, so its
triality dot is at most 6√2; the dense and single-row caps are 10√2
and 2√2. The new five-row layer adds 885,760 points without removals.
The pinned replay checks all pairs of tag patterns and the inner code.
