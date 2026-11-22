# GF(8) norm pole replacement in dimensions 63–67

| Dimension | Previous lower bound | New lower bound | Gain |
|---:|---:|---:|---:|
| 63 | 184,426,297 | **184,427,081** | 784 |
| 64 | 335,795,354 | **335,796,250** | 896 |
| 65 | 336,256,796 | **336,257,692** | 896 |
| 66 | 342,389,574 | **342,390,470** | 896 |
| 67 | 355,195,404 | **355,196,300** | 896 |

Replace each ordinary antipodal coordinate pole pair by sixteen signed
vectors. Each replacement has absolute value 6 on that axis, 2/3 in
all other 63 ordinary positions, and zero tail, so its squared norm is 64.
All previously certified sparse and lifted points are retained.

Index the ordinary coordinates by GF(64), with the zero last. For a
marked position a the field norm (x+a)^9 to GF(8) labels one position
by zero and nine by each nonzero label. Its binary coordinate functions
span a three-dimensional code inside the doubly-even BCH radical D.
Every nonzero word has weight 36 and vanishes at a. Translate the eight
words to an exterior coset of D in its extended BCH parent P, adding
their complements for the opposite heavy sign. The 64 translated groups
are disjoint, and avoid the incumbent exterior BCH cap E and all three
shifted Kerdock owner groups. The 64D dense stock can be realized as D
with the same cardinality and distance at least 16, regardless of the
inherited sparse signs. This is also the dense layer in the exterior-cap
construction.

Same-axis equal-heavy-sign pairs attain dot 32; opposite heavy signs
give at most −8. For different axes the distinct sign words are in P,
hence have distance at least 14. Even granting agreement on both heavy
positions and losing two of the guaranteed differences there, the dot
is at most 8+(4/9)(62−24)<32. Against a flat body in D or a lifted
flat body, distinctness and distance 14 give raw dot at most 88/3.
The lift scales this down and contributes no tail term. Against the old
signed 16-, 4-, and 9-site layers, triangle inequalities give
2(6+15·2/3), 4(6+3·2/3), and (8/3)(6+8·2/3), each at most 32.
They apply to the recent Netto substitutions too.

## The balanced 63D section

Relabel the GF(8) subfield octad to the pinned physical octad. Its
projection on D has rank seven and consists of all even words;
there are 70·2^21 dense sign words with exactly four minus signs, the
same count as the old section. Mark only coordinates outside that
octad. The norm restricted to the octad splits into four pairs, which
receive opposite signs under a suitable exterior translation.
Norm-word changes preserve this pairing balance. Thus the 896
spikes lie in the same hyperplane and replace exactly its 112 poles;
no other stock points or their coordinates are changed.

Run `node verify.cjs` in the repository package. The checker pins the
inherited claims and original section mask, and replays the finite
field, cosets, section balance, and all new pairs using exact arithmetic.
For this fixed one-heavy profile and heavy sign, the eight-word group
already attains the length-63 distance-36 Plotkin bound; this does not
cap other amplitude profiles or replacement geometries.

An independent `python verify.py` (with NumPy for the pinned BCH module)
also certifies the same five totals. It constructs different valid
translates, choosing a radical coset distinct from the shifted Kerdock
caps, then checks the physical balance and new pair products separately.
