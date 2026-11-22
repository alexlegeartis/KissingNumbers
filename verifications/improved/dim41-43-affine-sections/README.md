# Affine kissing configurations in dimensions 41–43

| Dimension | Parent | Complete families | Points in complete families | Individual points with a body | Pure head points | Total | Gain over previous repository bound |
|---:|:---|---:|---:|---:|---:|---:|---:|
| 41 | P48q | 151 | 1,460,064 | 0 | 2 | **1,460,066** | 30 |
| 42 | P48q | 184 | 1,991,469 | 0 | 6 | **1,991,475** | 30 |
| 43 | P48p | 119 | 2,917,905 | 3,165 | 3 | **2,921,073** | 1,270 |

## How the constructions were found

In 41D a search over the 729 affine centres of the complement of a scaled E7 in P48q produced a section with 1,164 points in its smallest body layer. In 42D the search kept the head coordinates and varied the body classes under E8 Weyl transformations. A difference class with an empty short shell permits an additional family. An orbit of 161,280 ordered pairs contains a choice that contributes 30 additional points. Both witnesses below use the rank-40 complement as a common coordinate system. [Earlier 41–43D constructions](../dim41-46-p48-heads/).

In 43D a discrete search moved the class pattern through the E8 Weyl group across the three P48 complements. The winning orientation in P48p puts 1,390,020 points in the twelve families of body squared radius 16, gaining 1,080 points in that layer. A joint integer selection moved complete families to rational heads and admitted a 61st family of 495 points. A second conflict graph reassigned individual small-body points and selected the pure head directions. The pinned upstream data identify the P48 sections used to start the search. [Source lattices](https://github.com/wcgbg/sphere-packing/tree/a0bc7f51444e42c3cbae0d2a910039cd50013500/packings_out/records/dim40-delta41).

## Coordinates and verification

Each archive stores a positive definite integral Gram matrix \(G\) of rank 40. A full family consists of every integer vector \(x\equiv c\pmod 3\) with \(x^TGx=3r\), paired with the rational head \(z/D\). Use the body metric \(G/3\) and normalize the resulting vector to unit length. Individual small-body points use the same metric, their own common head denominator, and \(r=10\). A pure head point is the unit vector in its specified direction. The three rational head spans have ranks one, two and three.

For two vectors with body radii \(r,s\), heads \(z/D,w/E\), and \(N=rD^2+z^Tz\), \(M=sE^2+w^Tw\), their normalized inner product is at most one half exactly when

\[
  2x^TGy\ \leq\ \left\lfloor
       \frac{\lfloor\sqrt{9NM}\rfloor-6z^Tw}{DE}
       \right\rfloor .
\]

The verifier enumerates each complete shell from the supplied Gram matrix, residue and radius. For two complete families it enumerates every possible nonzero difference shell that could violate the displayed inequality. Evenness of \(G\) and \(Gc\equiv0\pmod3\) fix each difference norm modulo 18. The verifier checks zero differences separately. For an individual point and a complete family it also tests whether adding an enumerated short difference to that individual point produces a member of the specified family. Direct integer products cover all 5,007,030 pairs of individual points in 43D and all comparisons with pure heads. Integer ellipsoid enclosures and exact leaf tests certify the enumeration. [Shared enumerator](../dim40-affine-centre-refill/enumerate_core.cpp).

From the repository root, with Python, NumPy and a C++17 compiler installed, run

```sh
python verifications/improved/dim41-43-affine-sections/verify.py --workers 12
```

Supply a dimension such as `43` to replay just that configuration. The faster `final.py` reads the population ledger. [Recorded exact replays](verification.recorded.json).
