# A lattice improvement in dimension 73

25 October 2026. The explicit lattice in this package proves

\[
\tau_{\mathrm{lattice}}(73)\ \ge\ 6\,218\,188\,144.
\]

This improves the **previous repository bound of 6,218,187,616 by 528**.
It improves the direct-sum lattice bound, 6,218,175,602, by 12,542.
The comparison is with the repository checkpoint on 25 October 2026.
A current comprehensive table for dimension 73 was not located; global
priority is not established by this package.

The lattice has minimum squared norm **8** and Gram determinant **5/2**.
`gram73_twice.json` is an explicit 73 by 73 integer matrix **twice** its
Gram matrix. Equivalently that integer matrix itself defines a scaled
lattice with minimum squared norm 16 and the same kissing number.

## The construction

Let \(L=\Gamma_{72}\), with minimum squared norm 8 and 6,218,175,600 minimal
vectors. These parent facts are from
[Nebe's theorem](https://arxiv.org/abs/1008.2862) and the
[Nebe–Sloane catalogue entry](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/neb72.html).

We exhibit \(u\in L\) such that the coset \(u+2L\) has minimum squared
norm 22 and contains at least 6,272 vectors of that norm. Define

\[
L_{73}=\{(x+k u/2,\ k\sqrt{5/2}):x\in L,\ k\in\mathbb Z\}.
\]

This is a full rank lattice. Its layers have the following minima:

- For \(k=0\), the nonzero minimum is 8, inherited from \(L\).
- For \(k=\pm1\), the body lies in \(u/2+L\), with minimum \(22/4=11/2\).
  Including the last coordinate gives \(11/2+5/2=8\).
- For \(|k|\ge2\), the last coordinate alone contributes at least 10.

Thus no lattice vector is shorter than the claimed contact vectors, and

\[
\tau(L_{73})\ge 6\,218\,175\,600+2(6\,272)
                 =6\,218\,188\,144.
\]

Writing \(G\) for the catalogue Gram matrix and \(u\) for a row of integer
coefficients, the Gram matrix of the displayed basis is

\[
G_{73}=\begin{pmatrix}
G&Gu^{\mathsf T}/2\\
uG/2&8
\end{pmatrix}.
\]

Its Schur complement is \(8-22/4=5/2\). Since \(\det G=1\), this proves
both positive definiteness and \(\det G_{73}=5/2\).

## What changed

The previous coset construction retained only vectors whose three Leech
block costs were at most 5. Completing the two missing types adds:

| Sorted block costs | Previous vectors | New vectors |
|---|---:|---:|
| (2, 2, 7) | 0 | 24 |
| (2, 3, 6) | 0 | 240 |
| (2, 4, 5) | 848 | 0 |
| (3, 3, 5) | 1,704 | 0 |
| (3, 4, 4) | 3,456 | 0 |
| Total | 6,008 | 264 |

Each additional coset vector gives one contact in each of the two layers.
The existing seed and its certified nearestness come from the repository's
`dim73-88-gamma72-cosets` package. The new work completes its omitted block
types and gives the explicit lattice interpretation and Gram matrix.
The two other tested seed cosets also yielded 6,272 vectors, from previous
subsets of 6,000 and 5,998; no equivalence claim is made for these cosets.

## Verification

From this directory, with Python and NumPy installed:

```text
python verify.py
python verify.py --rebuild
```

The first command takes about two seconds on the research machine. It
checks all 6,272 representatives, their norms, distinctness and common
coset; the exact matrix identities; the Gram matrix; and the absence of
shorter vectors using 58,608,516 pairs of block labels.

The second command also compiles the supplied C++ sources with `g++` and
reruns all three exhaustive block enumerations and the 292 targeted
completion enumerations. It takes about a minute. Both commands were run.
The interval bounds are checked against an exact rational LDL
decomposition. No time-limited search is accepted as a completed enumeration.

`certificate.npz` stores the integer matrices, 6,272 integer coefficient
rows, and the three complete block lists. `verification.json` records the
quick replay; `verification_rebuild.json` records the full enumeration replay.
`result.json` records the comparison and source facts.
`catalogue_gram72.json` was independently retrieved from the published
catalogue on 25 October; the verifier compares all 5,184 entries with the
parent Gram matrix used in the certificate.

The bound is for an actual lattice, not just an antipodal spherical code.
The 6,218,175,600 parent vectors are supplied by Nebe's theorem, rather
than materialized as a multibillion-row file. We claim a lower bound on
the new lattice's kissing number; this package need not establish that
the displayed coset list is exhaustive to prove the improvement.

## Minimum certificate in more detail

The supplied matrices \(F,A,B,T\) satisfy \(A+B=I\), \(AB=2I\),
\(AFA^T=2F\), and

\[
T=\begin{pmatrix}I&I&B\\0&A&A\\0&0&2I\end{pmatrix},\qquad
T\operatorname{diag}(F,F,F)T^T=2G.
\]

Here \(F\) is the Leech Gram matrix. A coefficient row \(x\) has block
coordinates \(w=xT\), and its norm is half the sum of the block norms.
In the selected parity class write \(w_i=q_i+2z_i\), with the three fixed
\(q_i\) supplied by the first two rows of \(A\). Define
\(c_i=w_iFw_i^T/4\). For the parity class in this certificate all these
costs are integers, and the complete enumerations give \(c_i\ge2\).

Every vector in the parity class necessarily satisfies, modulo 2,

\[
z_0B=z_1B,\qquad z_2=h+z_0A+z_1,
\]

where the verifier derives \(h=12,625,758\) from the stored centre.
Because \(L\) is even, all norms in \(u+2L\) are 22 modulo 4. A shorter
vector would therefore have \(c_0+c_1+c_2\le9\). Since each cost is at
least 2, each block then has cost at most 5. The three complete block
lists suffice to rule this out. The verifier evaluates all the necessary
label matches and finds zero shorter triples.

The block lists are regenerated using outward rounded interval arithmetic;
their LDL intervals enclose exact rational values. Accepted leaves are
tested using integer quadratic forms. The explicit new heads are checked
directly using \(G\), independently of their discovery procedure.
