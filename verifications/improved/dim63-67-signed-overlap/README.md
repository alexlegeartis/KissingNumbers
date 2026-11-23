# Signed overlap-10 construction in dimension 63–67

Let C be the 80 cosets of binary RM(2,6) in the literal cubic palette in
`cubic80.txt`. The exact minimum distance is 16. Take all sign words
from C (335,544,320 vectors), the 166,656 signed supports on affine
two-flats, and coordinate poles. The previous construction included
37,996 weight-16 support classes (plus the stated shelves with 2 tail
coordinates). It used only support pairs of intersection at most eight.

Here we add 2,354 distinct weight-16, zero-tail supports. Their literal
coordinates and 30-bit-sized check records reside in `signed_witness.txt.gz`.
The new supports arose by searching five cubic translates, but the verifier
needs no claim concerning their provenance: it checks every geometric
pair passing the overlap-eight automatic threshold.

For each support S, after fixing a coordinate order, take its sign words
from an affine binary [16,11,4] extended Hamming code. The proof witness
lets this code vary with S. A length-16 binary check space which contains
the all-ones row and whose other nonzero rows have weight eight extends to
the dual of an extended Hamming code: the row-pattern classes of a basis
have equally sized fibres (by binary Fourier inversion). Split each class
in half successively, obtaining four independent, balanced coordinate
bits. Extend the check syndrome arbitrarily on the new basis; the sign
coset has 2048 members and min distance at least four because the five
check bits distinguish all coordinate columns and include overall parity.

The verifier confirms that every pair of supports meeting in 9 or 10
positions has a **common weight-eight subset with opposite sign parity**.
Their sign vectors therefore disagree in at least one of the shared
positions, giving dot ≤ 10−2 = 8. All other pairs of supports immediately
give dot ≤ 8 and hence cosine ≤ 1/2. All sparse norms squared are 16;
within one coset, Hamming distance four gives dot ≤8. The dense and
affine-flats layers have norms squared 64 and 4; worst mixed dots are
16 (dense/sparse), 4 (dense/affine or sparse/affine). The coordinate poles
are also automatically compatible. Every support uses norm 16 even when
two coordinates are in the tail.

For dimension 63 intersect the 64D code with a coordinate-equality plane.
The retained fraction on supports containing both coordinates is 1/2
because no two-column parity is in the five-check dual. Supports containing
neither are retained in full. If both are present in a dense sign word,
RM(1,6) makes the pair parity balanced; 124 pure poles survive.

One-command exact replay:
```sh
g++ -O3 -march=native -std=c++20 verify.cpp -o verify && ./verify
```

Output rows exclude any polar 65 witnesses, which may be added separately
because all new supports are zero-tail and compatibility follows by the
body ℓ¹ bound. The 68D row from this same verifier is not a new record.

The portable entry point is `python3 verify.py` on an x86-64 machine
with BMI2 and a C++20 compiler. It unpacks the pinned witness. In 65D
we retain the three BCH polar cells on each side from the
[separate polar proof](../dim63-66-internal-tags-bch/), and replay those
cells again here. Their 393,216 points are valid against the new
equator by its sparse L1 bound, with no need to keep the earlier
internal amplitude-four tags. The published counts are:

| Dimension | Lower bound |
| --- | ---: |
| 63 | 216,972,548 |
| 64 | 418,347,904 |
| 65 | 418,741,122 |
| 66 | 426,343,300 |
| 67 | 442,334,086 |
