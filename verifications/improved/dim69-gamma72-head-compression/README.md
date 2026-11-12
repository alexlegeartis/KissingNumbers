# A 69-dimensional kissing configuration with at least 743,522,304 points

**Retained first-stage construction.** The [extended D7 moment package](../dim68-69-gamma72-extended-moments/) strengthens this geometry to **746,325,051** points, including two pure head poles. The original certificate and its historical gain below are preserved.

**13 November 2026.** This package proves

\[
\boxed{\tau(69)\ge743\,522\,304.}
\]

The preceding repository certificate gives 689,452,922. The increase is
**54,069,382 points, or 7.8424%**. This is an unrestricted spherical construction
obtained by deforming slices of a lattice shell. The population proof establishes
that a favorable orientation exists; it does not list all resulting coordinates.

## Reproduce

Python 3 and its standard library suffice. From this directory:

```sh
python verify.py
python verify_orthonormal.py
python rejection_controls.py
python final.py
```

The two verifiers reconstruct the proof using different coordinates and moment
formulas. The first uses the integral Gram profiles and a Gaussian recurrence.
The second uses orthonormal coordinates, signed permutations and triality, and
the closed formula for spherical monomial moments. Neither imports an optimizer,
NumPy, SymPy, or discovery code. `final.py` supplies the registry through
`results()` only after the independent exact verification succeeds.

## 1. Published starting lattice and an actual normal sublattice

Let \(L=\Gamma_{72}\) be Nebe's extremal even unimodular lattice. We use these
published facts: its squared minimum is 8, its minimal shell \(X\) has
6,218,175,600 vectors, and that shell is a spherical 11-design.

The complete catalogue Gram matrix and four integer coefficient vectors are
included in `embedding.json`. Direct integer matrix multiplication checks that
the selected normal vectors have Gram matrix

\[
G=\begin{pmatrix}
8&-4&-4&-4\\
-4&8&0&0\\
-4&0&8&0\\
-4&0&0&8
\end{pmatrix}.
\]

Thus the normal span \(W\) is an actual copy of the scaled root lattice
\(2D_4\) inside \(L\). No unproved embedding or hypothetical extremal lattice
is used. Put \(V=W^\perp\), so \(\dim V=68\).

For clarity, the 11-design fact follows from the standard extremal-lattice
theorem. For a homogeneous harmonic polynomial of positive even degree
\(d\le10\), the weighted theta series has weight \(36+d\) and vanishes to
order at least four at the cusp. Dividing by \(\Delta^4\) would give a
holomorphic modular form of weight \(d-12<0\), so the series is zero.
Antipodality treats odd degrees. Theorem 4.1 in Nebe's survey below gives this
argument in general. The verifiers take these published lattice and design
facts as inputs; they do not reprove the modular-form theorem.

## 2. Replace four normal coordinates by one adjustable coordinate

Write a shell vector as \(v=b+p\), with \(b\in V\), \(p\in W\). Choose
orthonormal coordinates on \(W\), and write \(p=h/4\). For each selected
profile \(h\), choose a rational height \(t_h\), and define

\[
q_h=\frac{\lVert h\rVert^2}{16},\qquad
\alpha_h^2=\frac{8-t_h^2}{8-q_h},\qquad
F(v)=(\alpha_h b,t_h)\in V\oplus\mathbb R.
\]

All image vectors have squared norm 8. The selected ten profiles and heights
are as follows; the coordinate order matches `orthonormal_certificate.json`.

| Profile \(h\) | \(q_h\) | Height \(t_h\) |
|---|---:|---:|
| \((-1,1,1,1)\) | \(1/4\) | \(1/8\) |
| \((-2,-2,-2,2)\) | \(1\) | \(-1\) |
| \((-1,-1,1,1)\) | \(1/4\) | \(-2/15\) |
| \((0,0,0,2)\) | \(1/4\) | \(0\) |
| \((1,1,1,1)\) | \(1/4\) | \(1/2\) |
| \((-1,-1,-1,1)\) | \(1/4\) | \(-1/2\) |
| \((-2,0,0,2)\) | \(1/2\) | \(-4/15\) |
| \((0,0,2,2)\) | \(1/2\) | \(4/15\) |
| \((1,1,1,3)\) | \(3/4\) | \(1/2\) |
| \((2,2,2,2)\) | \(1\) | \(1\) |

This replaces a four-dimensional projected component by one coordinate while
allowing different slices to receive different heights and rescalings. The
choice allows five populous profiles of squared projection norm \(1/4\) to
coexist with five further slices.

## 3. Every angle inequality is exact

Distinct minima \(v,w\in X\) satisfy \(v\cdot w\le4\), since
\(\lVert v-w\rVert^2\ge8\). For profiles \(h,k\), put
\(q_{hk}=h\cdot k/16\). Their transformed inner product is at most

\[
\alpha_h\alpha_k(4-q_{hk})+t_ht_k.
\]

Every pair of selected profiles satisfies

\[
\alpha_h\alpha_k(4-q_{hk})+t_ht_k\le4.
\]

The verifiers first check the signs of both sides and then prove the equivalent
rational inequality

\[
(4-t_ht_k)^2-\alpha_h^2\alpha_k^2(4-q_{hk})^2\ge0.
\]

This checks **all 55 same-slice and cross-slice pair types**. Six are tight; the
smallest positive squared margin is \(1/180\). Therefore the normalized image
is a kissing configuration. This also proves injectivity: coincident norm-eight
images would have inner product 8, contradicting the bound 4 for distinct
source vectors.

## 4. Complete profile domain and valid orientation averaging

Integral products with all roots of \(2D_4\) imply that the coordinates of
\(h=4p\) are either all even or all odd. Apart from the 24 known normal roots,
the lattice minimum gives

\[
|h_i+h_j|\le8,\qquad |h_i-h_j|\le8\quad(i\ne j).
\]

Equivalently, the two largest absolute coordinates sum to at most 8.
Exactly **1,681** profiles satisfy these conditions. The known roots have
profiles \(\pm8e_i\pm8e_j\), and each contributes exactly one source vector.

The root lattice admits signed coordinate permutations and the orthogonal
Hadamard map

\[
\frac12\begin{pmatrix}
1&1&1&1\\1&1&-1&-1\\1&-1&1&-1\\1&-1&-1&1
\end{pmatrix}.
\]

They preserve the root lattice and hence take the chosen normals to actual
integer combinations of lattice vectors. Their generated group is finite
because it acts faithfully on a finite spanning set of roots. Average the
actual profile distributions over this group. The resulting counts are constant
on **16** profile orbits and obey all design moments.

This argument does **not** require the ambient lattice or its shell to be
invariant under these transformations. Each transformed normal tuple exists in
the same lattice, and the average of their valid distributions is well defined.
The ten-slice population averaged over orientations is at most the largest
actual ten-slice population. Consequently a lower bound on the average is an
existence bound for one actual orientation.

## 5. Rational population certificate

For nonnegative integers \(a_1,\ldots,a_4\), with \(s=\sum a_i\le5\), the
11-design supplies

\[
\sum_{v\in X}\prod_{i=1}^4 h_i(v)^{2a_i}
=6\,218\,175\,600\,
\frac{128^s\prod_i(2a_i-1)!!}
{72\cdot74\cdots(72+2s-2)}.
\]

Empty products and \((-1)!!\) are one. Subtract the 24 known roots, collect
the 16 orbit variables, and clear row denominators. Ten independent rows suffice
for this certificate. Write the resulting equations as \(Ax=b\), with
\(x\ge0\). The ten-slice objective is \(c^Tx\).

The orthonormal checker reconstructs a rational vector \(y\) with denominator
1,451,520 and numerator vector

```text
0, 12791808, -6743552, 1314096, 82752,
-29934, -160626, -3916, -3333, 10305
```

It checks \(A^Ty\le c\) entry by entry and obtains, exactly,

\[
b^Ty=743\,522\,304.
\]

Thus every feasible averaged distribution has ten-slice population at least
that value. A favorable actual orientation has at least as many points, and
Sections 2–3 turn those points into a valid code in dimension 69. There is no
rounding loss and no use of an optimizer's numerical status as a proof.

## Files and proof scope

- `certificate.json`, `verify.py`: exact Gram-coordinate certificate and checker.
- `orthonormal_certificate.json`, `verify_orthonormal.py`: independently rebuilt
  orthonormal certificate and checker.
- `embedding.json`: the pinned catalogue matrix and actual lattice normals.
- `final.py`: exact verification and current-record driver.
- `source_review.json`: dated baseline, primary sources and comparison scope.
- `verification.json`, `verification_orthonormal.json`: exact replay receipts.
- `rejection_controls.py`, `rejection_controls.json`: package-local semantic
  corruption tests for both checkers; two valid baselines and 26 rejections.
- `INDEPENDENT_REVIEW.md`, `independent_review.json`: third adversarial
  mathematical and arithmetic review, including the corrected `-O` guard.
- `SHA256.json`: content hashes of the final proof package.

Discovery scripts and the bounded neighboring-dimension searches are preserved
in `../../../research/dim32plus_20261113/`. This proof applies to the full selected
families using moment counting; no assertion of an explicit favorable orientation,
full coordinate enumeration, or global optimality is made.

## Primary references

1. G. Nebe, *An even unimodular 72-dimensional lattice of minimum 8*,
   J. reine angew. Math. 673 (2012), 237–247.
   <https://arxiv.org/abs/1008.2862v3>.
2. G. Nebe and N. J. A. Sloane, *The Lattice Gamma72*, Catalogue of Lattices:
   <https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/neb72.html>.
3. G. Nebe, *Boris Venkov's Theory of Lattices and Spherical Designs*,
   Theorem 4.1:
   <https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/venklat.pdf>.
4. Preceding 69D proof in this repository:
   `../dim68-69-gamma72-cross-sections/weyl_averaged/`.
