# Post–Reed–Solomon partial fibres in 94–96 dimensions

| Dimension | Previous bound | Gain | New bound |
|---:|---:|---:|---:|
| 94 | 74,142,076,620 | 121,664 | **74,142,198,284** |
| 95 | 117,091,749,956 | 121,664 | **117,091,871,620** |
| 96 | 209,293,999,424 | 33,056 | **209,294,032,480** |

In 94D and 95D, append 148 six-row triality modes with 3,802 explicitly
listed axis words between them. Each takes all 32 even sign choices.
All preceding vectors are retained, including the 164-word codes of the
Reed–Solomon repacking. The new supports are entirely in the common
88-coordinate body, so the same geometry works in both dimensions. The
ordinary/paired BCH, D4 matching, flat RS, grouped root, amplitude-four RS,
dense and shifted Golay interfaces use the published purely geometric
bounds, which are independent of the axis labels.

In 96D, instead replace the axis words in 76 modes that still use 128
words in the prior bound. Each replacement is strictly larger and gains
a net 1,033 axis words in total. These are deliberately different modes
from the 1,334 that were already repacked. Every changed mode has at least
one active partial-fibre neighbor, but none has any active full-mode
neighbor. This makes the replacement decisions independent of one another,
and the verifier checks the exact projections against the literal 60
retained partial fibres. No other old mode or geometric layer is changed.

When two modes share q physical rows and have the same frame on r of
those rows, their raw product is at most 2(q+r)−4d, where d is the
Hamming distance after projecting both words to the common matched
frames. This is an integer upper bound because a distinct axis of one
24-cell frame is orthogonal and different frames cost at most 2.
The cap is 12; the verifier checks every active interface using direct
projected axis-word sets, without solver assumptions. Each new mode has
minimum axis-word distance 3 internally, and distinct even sign choices
have distance at least 2, also giving the cap 12.

Run `python3 verify.py` with assertions enabled from the repository.
The sibling files are pinned by Git blob hashes. Every line of
`witness94_95.txt` lists the three disjoint frame masks, count and axis
words. Every line of `witness96.txt` gives a zero-based predecessor
mode index, count and axis words. Digits are stored in two-bit fields in
increasing physical row order. The discovery scripts selected
intersections of 164-word or 128-word codes and optimized their union;
these scripts and the numerical search are not needed for verification.
