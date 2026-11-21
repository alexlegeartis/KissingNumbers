# Two asymmetric latitude pairs in dimension 89

**15,125,338,114 kissing points**, up **1,073,741,824** from the preceding repository value 14,051,596,290 (about **7.64%**). This is a lower bound from explicit fibers and the previously verified sparse layers, not an assertion of optimality.

## Construction and proof

Use the [exact 88D and 89D Z8 parent](https://github.com/alexlegeartis/KissingNumbers/tree/27247d2c57e82d589f688b2996ec2dfabf8ac4fa/verifications/improved/dim88-89-z8-prefix-anticode). The old mixed-latitude configuration expands only one of the three crosspolytope colors, with 5+4+4 fibers. The new construction expands **two colors simultaneously** by swapping their large and small signed heights. Delete the first two four-dimensional blocks as before; let `b` denote the second deleted symbol, with first deleted symbol 0 throughout. All body vectors have squared norm 88.

| Color | `b` at height 0 | `b=3`, positive | `b=7`, negative |
|---|---|---|---|
| 0 | 0,1,2 | √8 | −√8/4 |
| 1 | 0,1,2 | √8/4 | −√8 |
| 2 | 0,1,2,3 | (0), already counted | absent |

Normalize each vector by its **individual** radius (sqrt{88+h^2}). By the pinned source minimum, body products are bounded by 40 within the same fiber, by 44 across distinct nonopposite fibers in a color, by 48 in the two opposite pairs, and by 44 across colors. Every squared radius is at least 88. Within a fiber (40+h^2le(88+h^2)/2) since (h^2le8). The two same-color opposite pairs have raw total dot at most (48-2=46). Cross-color nonzero-height pairs of the same sign also have raw total dot at most (44+2=46); the respective radii in both cases are (sqrt{96}) and (sqrt{177/2}), and (4cdot46^2=8464<8496=96cdot177/2). Other distinct fiber pairs have raw product at most 44. The two pure-head poles remain compatible: maximum (h^2=8) and (4h^2le88+h^2). Strict pair cap below 1 ensures no directional duplicates.

Each fiber retains exactly (8^{10}) words because the two-deleted-coordinate map is onto and the triangular odd pivots give injection; this is checked on the pinned exact generator. Their total with poles is (14cdot8^{10}+2=15,032,385,538). The previously published weight-eleven and weight-ten Golay sparse layers are unchanged. Their heads vanish, and each occupied row gives body dot at most 4 against all three dense colors, exactly as before. For (w=11) or (10), (4(4w)^2le(8w)(88+h^2)), so they are compatible with every new dense fiber. They add (629cdot1024cdot144+11cdot256cdot72=92,952,576) points. Distinct support and squared body radii exclude duplicates with the dense fibers or poles.

## Replay

Run `python3 verify.py` (also requires Node.js for the unchanged sparse proof). It independently reconstructs the three alphabets and all 91 inter-fiber pair bounds, checks same-fiber distances, poles, sparse–dense coupling, SHA-bound source generator, and arithmetic. The inherited published Z8 homogeneous-minimum-24 proof remains an explicit premise, not inferred from rank. `abbreviated_parent/` carries the original checker and data for the unchanged sparse layer. The new dense witness replaces, rather than adds to, the previous dense layer.

`claims.json` stores the new total, `verify.log` the exact replay output, and `z8_generator.json` the pinned source generator.
