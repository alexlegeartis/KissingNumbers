# Mixed three-tetrad fibres in 94D and 95D

| Dimension | Previous lower bound | Added points | New lower bound |
|---:|---:|---:|---:|
| 94 | 74,147,177,220 | 279,552 | **74,147,456,772** |
| 95 | 117,096,850,556 | 279,552 | **117,097,130,108** |

All preceding points are retained. The new shell has six modes with
1,152 nonempty axis fibres. Every support avoids rows zero and one, so
the same body fibres extend by zero in the head in both dimensions.

Each mode comes from two binary Golay octads with intersection four.
The three tetrads of their union, in the order (a\setminus b),
(b\setminus a), (a\cap b), use the three D4 matching frames. Sort
the physical rows within each tetrad. Every occupied row has two
entries of magnitude 2 and squared norm 8; the declaration of axes
is a permuted 4,096-word dodecacode of quaternary distance six, with
one XOR shift. Distinct axis words then give at most six aligned rows.
For the same axis, the signs have minimum binary distance three.

The 1,064 actual mask-3 shifted predecessor modes are loaded by blob
hash. All but fourteen geometric interfaces are safe before phases;
on those fourteen, q occupied rows are common and r have the same
frame, with q+r=16. The only potentially dangerous old axis words
differ in zero or one of the r matched slots. Their nonzero row products,
in units of four, are signed 2 on equal axes and signed 1 on differing
frames. At raw absolute sum 14, all active signs must align to exceed
12; at sum 16, they must either all align or differ on a single
different-frame row. Remove all such achievable patterns using the
old shortened Hamming syndrome reduced modulo the erased columns.
The [96D replay](../dim96-mixed-tetrad-axis-sign-fibres/) uses the same
exact projection rule.

For each remaining sign fibre, compare the best standard Hamming
coset to an adapted XOR syndrome. The latter assigns colors to the
twelve sign-coordinate labels and a zero sentinel: join a pair when
the safe fibre contains a difference in those two coordinates, or
join to the sentinel if the coordinate can flip alone. A deterministic
saturation-degree coloring gives labels such that a fixed XOR syndrome
excludes every one- and two-bit difference. Retain the largest coset
from the better recipe. This calculation does not assume that a safe
sign fibre is affine.

Every new-new interface and every interface with the 328 current
coordinate-aligned monochromatic modes is safe by q+r≤12. The 85 short
matching modes of weight w satisfy (q+r)^2≤12w, the exact squared cap
for their normalized mixed products. The dense, RS, BCH, grouped root
and smaller six-row triality stock is protected by the same per-row
geometry as the previous shifted extension.

Run `python3 verify.py` with assertions enabled and g++ available.
The wrapper compiles the C++17 integer replay in a temporary directory.
Each line of `witness.txt` gives two octad masks, a permutation of
the 12 slot digits 0–9ab, a 24-bit axis XOR shift, and the verified
point count. No external packages, processor-specific compiler flags,
or optimizer are needed.
