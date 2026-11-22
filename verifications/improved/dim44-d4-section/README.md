# A 4,940,961-point kissing configuration in 44 dimensions

Choose four norm-six vectors generating a scaled copy of the D4 root lattice in the extremal even unimodular lattice P48m. Project seventeen affine shells to their 44-dimensional orthogonal complement and normalize each vector. Varying the four normals made three projected shells of squared radius three empty. This gains **1,680 points** over the previous repository bound of 4,939,281. The parent Gram matrix was obtained by an integral change of basis from the pinned P48m certificate. The replay verifies extremality from the supplied matrix. [Source Gram matrix](https://github.com/wcgbg/sphere-packing/blob/a0bc7f51444e42c3cbae0d2a910039cd50013500/packings_out/records/dim45-kissing7379838/P48m-d914e7d8.cert.txt).

| Projected squared radius | Number of families | Points per family | Subtotal |
|---:|---:|---:|---:|
| 6 | 1 | 4,127,904 | 4,127,904 |
| 5 | 9 | 90,156 | 811,404 |
| 4 | 1 | 825 | 825 |
| 11/3 | 3 | 142 | 426 |
| 11/3 | 3 | 134 | 402 |
| **Total** | **17** | | **4,940,961** |

## Reconstructing the vectors

The rows of `normals` in [ambient.npz](ambient.npz) are lattice vectors in the integral basis with Gram matrix `gram`. Write their matrix as \(B\), their Gram matrix as \(K\) and its determinant as \(D=324\). For a parent vector \(x\) with integer inner products \(v=xGB^T\), its projected vector is \(x-vK^{-1}B\). Its class label is \(v\operatorname{adj}(K)\bmod D\).

[The family ledger](families.npz) gives the seventeen class labels, squared radii in units of \(1/D\), and populations. For each label choose any representative \(v\) minimizing \(vK^{-1}v^T\), and write the minimum as \(a\). Project every lattice vector with those four inner products and parent squared norm \(r+a\), where \(r\) is the prescribed projected squared radius. Representatives differing by a vector of the normal lattice give identical projected sets.

For different projected points in the same family the squared distance is at least six. For different families, lift a nonzero difference using a shortest normal representative of its difference class. If that representative has squared length \(a'\), the projected squared distance is at least \(6-a'\). The verifier checks \(r+s-(6-a')\leq\sqrt{rs}\) by integer arithmetic for every cross-family pair. The two families with the same label have different squared radii.

## Why the populations are exact

There are 324 classes, or 164 up to sign. Every class has a normal representative of squared length at most three. This also proves the normal lattice is primitive. Any additional parent vector in its real span could be translated to a nonzero vector of squared length at most three, contradicting the verified minimum six. The spherical 11-design equations on the parent norm-six and norm-eight shells give 160 independent ordinary moment equations per shell. Sixteen chosen families come from shell six and one comes from shell eight. A degree-twelve harmonic theta series has weight 36 and vanishes through degree two, hence is a multiple of the modular discriminant cubed. Its degree-four coefficient is −72 times its degree-three coefficient. Four such relations determine the remaining freedom on shell eight once the shell-six populations are fixed. [Extremal lattice shells and spherical designs](https://www.mathnet.ru/php/archive.phtml?jrnid=mmj&option_lang=eng&paperid=245&wshow=paper). [Classical spherical theta series, Theorem 5.18](https://publications.rwth-aachen.de/record/685470/files/685470.pdf).

[moments.py](moments.py) establishes the ranks by elimination modulo 1,000,003 and checks the supplied affine formulas over the integers. Four explicitly counted shell-six profiles complete the system. Let \(a,b\) be the first normal and the negative of the second. Both have squared norm six and \((a,b)=3\). Spherical moments give exactly 1,656 shell-six vectors with inner products \((x,a)=3\) and \((x,b)=0\). The archive supplies that many distinct vectors satisfying both equalities. Their inner products with the other two normals give the four needed counts, **0, 0, 0 and 134**.

For this choice of seventeen classes the size simplifies to \(4{,}940{,}961-70(c_1+c_2+c_3)\), where the first three counts are the populations of the squared-radius-three shells. This identity guided the discrete search over triples of mutually orthogonal neighbors of a fixed minimal vector.

With the first two normals fixed, [search.py](search.py) exhausts the remaining orthogonal pairs in the supplied 1,656-point slice. It tests 1,866 pairs and finds 111 attaining 4,940,961. The stored normals are one of those completions.

## Replay

From the repository root run

```sh
taskset -c 0-126 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verifications/improved/dim44-d4-section/verify.py --workers 12
```

[final.py](final.py) sums the seventeen verified populations for the repository ledger. The full replay uses Python, NumPy, SymPy and a C++17 compiler. It proves that the ambient Gram matrix is positive definite, even and unimodular, exhausts its vectors of squared norm at most four, checks the coordinate slice, and verifies all counts and angles. The [recorded output](verification.recorded.json) includes the four independently recovered calibrations.
