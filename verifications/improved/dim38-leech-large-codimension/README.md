# A large new lower bound for the kissing number in dimension 38

$$\tau_{38} \;\ge\; 591612$$

against the published record **566652** (H. Cohn, *Table of kissing number bounds*) — an
improvement of **+24960**, a factor 1.044.

(Earlier in this project the same dimension was improved to 570236 by feeding the 2026
constant-weight codes into the Edel–Rains–Sloane construction; see `../dim38-39/`. The
configuration here is a completely different one and beats that by +21376.)

```
python scripts/verify.py
```

runs in about ten seconds, uses only numpy, and decides everything in exact integer
arithmetic: 34 checks, none of them with a floating-point tolerance.

*The construction idea here — a class of Leech lines carried on a zero-sum triple of cap
directions — is the same one that gives dimension 27 and dimensions 73–95. What is new in
this package is running it at LARGE codimension, where the binding constraint flips from
"how large is one class" to "how much of the shell can be partitioned into classes".*

## The construction

This is the Leech-lattice "equator + caps + axis" construction — the one that holds every
record in dimensions 25 through 31 — run for the first time at a *large* codimension.
Split $\mathbb{R}^{38} = \mathbb{R}^{24}\oplus\mathbb{R}^{14}$ and take

| family | points | count |
|---|---|---|
| equator | $(y,0)$, $y$ a Leech minimal vector not used as a head | $196560-2c$ |
| caps | $(\sqrt{2/3}\,\hat u,\ \sqrt{1/3}\,\hat w)$ | $6c$ |
| axis | $(0,\hat a)$ | $A$ |

Here $W\subset S^{13}$ is a 60° spherical code — the 1932-point kissing configuration in
dimension 14 — partitioned into **zero-sum triples** $\{w_1,w_2,w_3\}$, $w_1+w_2+w_3=0$,
so that each triple's three directions are pairwise at $\cos = -1/2$. Each triple carries
a **class**: a set of Leech minimal vectors no two of which are at $\cos = 1/2$. Classes
on different triples are disjoint, and every vector used as a head is deleted from the
equator. Writing $c$ for the number of *lines* (antipodal pairs) covered by the classes,

$$N \;=\; (196560-2c) \;+\; 6c \;+\; A \;=\; 196560 + 4c + A .$$

All of this is forced by the 60° condition; the derivation is reproduced inside
`scripts/verify.py`, and each of the six kinds of pair reduces to a single integer
inequality.

## Why this is new

For $k = \dim - 24 \le 7$ — the range everyone has used — the binding constraint is the
number of triples: $\lfloor \tau_k/3\rfloor$ is 0, 2, 4, 8, 13, 24, 42, so at most
$42\times 496 = 20832$ vectors can ever be used as heads, and the classes are taken to be
maximum ones (496 vectors = 248 lines). The bound is then
$196560 + 4\,(\#\text{triples})\times 248 + \tau_k$, which reproduces the published
dimension-25-to-31 records exactly.

At $k = 14$ the constraint flips. There are up to $1932/3 = 644$ triples, far more than the
$\lceil 98280/248\rceil = 397$ maximum classes could ever fill, so the question becomes not
*how large is a class* but **how much of the Leech shell can be partitioned into classes at
all** — a graph-colouring question that had not been asked. The answer found here is:

> **all of it.** The 98280 lines of Leech minimal vectors partition into **644 pairwise
> disjoint classes**, none of which contains two lines at $\cos = 1/2$.

A **perfect** partition of the 1932 directions into 644 zero-sum triples was also found, so
the equator is *empty* and the cap layer alone has $6 \cdot 98280 = 589680$ points — exactly
the ceiling $3\tau_{24}$ of this family.

On top of that the axis layer is **not** empty, even though every direction is already used
as a cap direction. An axis point $(0,\hat a)$ only needs $\langle\hat a,\hat w\rangle \le
\sqrt3/2$ — that is 30° — whereas the covering radius of the 1932-point code turns out to be
**54.7°**. A rotated copy of the same 60° code, turned so that **not one** of its 1932
points falls inside a 30° cap, contributes all of them:

$$589680 + 1932 \;=\; 591612 .$$

The rotation is $Q = N/D$ with $N$ integral and $NN^{\mathsf T} = D^2 I$, so that check is
exact integer arithmetic too; how it was found is
[`scripts/axis_rotation.py`](scripts/axis_rotation.py) and the section below.

**What is left: nothing, inside this mechanism.** Both layers are at their ceilings now.
The cap layer cannot move: $6\cdot98280$ is $3\tau_{24}$ with the equator empty. The axis
layer is a $60^\circ$ code of $\mathbb R^{14}$, so it has at most $\tau(14)$ points, and it
has exactly $\tau(14)=1932$. So

$$591\,612 \;=\; 3\tau_{24} + \tau(14),$$

and the construction has no slack left at all. The one honest caveat is that $\tau(14)=1932$
is the best *known* kissing number in dimension 14 rather than a proven one: the ceiling is
$3\tau_{24}+\tau(14)$ with the true value, so anyone who improves dimension 14 improves this
by the same amount, and nothing else does.

### The 24 points that used to be missing

This read **1908** until 2026-08-27, and this README argued that the last 24 were out of
reach: the $30^\circ$ caps cover about $1.24\%$ of $S^{13}$, so a random rotation loses about
24 points and "best-of-$N$ would recover only a handful". Neither half survives. $1.24\%$ is
$24/1932$ — the loss of the one rotation in hand, read back as if it were a measure — and the
real figure is $2.89\%$: a $30^\circ$ cap on $S^{13}$ has relative measure $1.4978\times10^{-5}$
and $1932$ of them give $2.8937\%$, whence an expected loss of $55.9$, which is what 200 random
rotations actually show (56, sd 11); the independent-points model gives $e^{-56}$, not $e^{-24}$.
The second half was **a mechanism offered in place of a measurement**, and it was wrong.

$\mathrm{SO}(14)$ has 91 dimensions and the objective is differentiable once the maximum over
cap directions is softened, so it can be *descended* rather than sampled:

$$S_{ij}=\langle Qw_i,w_j\rangle,\qquad
  f(Q)=\sum_i \mathrm{softplus}\big(\mathrm{softmax}_j S_{ij}-4\sqrt3\big),$$

with $\partial S_{ij}/\partial Q = w_jw_i^{\mathsf T}$, so the gradient is a $14\times14$
matrix and the step $Q\mapsto \exp(-\eta\,\mathrm{skew}(GQ^{\mathsf T}))\,Q$ stays exactly on
the group. Started from the old rotation it reaches all 1932 in a couple of hundred steps.

The result then has to be made **rational**, because the check is integer arithmetic. Cayley
does that: $(I-S)(I+S)^{-1}$ is rational and orthogonal for every rational skew-symmetric
$S$, so the real solution gives $S=(I-Q)(I+Q)^{-1}$, which is rounded to
$\tfrac1{15}\mathbb Z$ and rebuilt exactly in fractions. The descent is run against a
**margin** of $0.15$ rather than against the threshold itself, because rounding moves every
inner product a little and a solution sitting on the boundary would not survive it; the
rounded rotation clears the threshold by $0.113$.

One line of `verify.py` had to grow with it. The old test squared the inner product and
compared with $48D^2$, which forced $48D^2<2^{62}$ and capped $D$ near $3\times10^8$; the new
$D$ has 18 digits. Squaring is unnecessary — $48D^2$ is not a perfect square, so $4\sqrt3D$
is irrational, and for an **integer** inner product "$\le 4\sqrt3D$" is exactly
"$\le\operatorname{isqrt}(48D^2)$". That is one big-integer square root, computed once, after
which the $1932\times1932$ comparison stays in `int64` and needs only $8D<2^{63}$.

## How the 644-class partition was found

It exploits a structural fact about maximum classes discovered in this project: every known
496-vector class consists of 124 orthogonal **pairs** whose differences span a 4-dimensional
totally isotropic subspace $R$ of $\Lambda/2\Lambda$. Consequently the whole class lies
inside

$$X_R \;=\; \{\ \ell \;:\; b(\ell,R) = 0\ \},\qquad b(x,y) = x\cdot y \bmod 2,$$

which is just **6120** of the 98280 lines — a 16-fold reduction. Local search restricted to
such an arena returns classes of 230–244 lines; the identical search on the full graph
returns 152. Sweeping over the 4096 arenas obtained from $R$ by the Golay sign-change
automorphisms of the Leech lattice, and falling back to unrestricted search once the arenas
are exhausted, covers everything in 644 classes.

## What the script verifies

Every step is exact; the only floating point in the file is in the printed commentary.

1. **The Golay code.** The twelve stored generators are expanded to 4096 words and the
   weight enumerator is checked to be $1+759x^8+2576x^{12}+759x^{16}+x^{24}$.
2. **The Leech lattice.** The 196560 minimal vectors are rebuilt from that code, checked to
   be distinct, of squared norm 32, and to have the Leech inner-product distribution
   $1,4600,47104,93150,47104,4600,1$ from a fixed vector.
3. **The directions.** 1932 distinct vectors of squared norm 8 whose maximum off-diagonal
   inner product is 4 — i.e. a genuine 60° spherical code (this is Cohn's published
   dimension-14 kissing configuration).
4. **The triples.** All 644 are pairwise disjoint, each sums to zero, and each has all three
   pairwise inner products equal to $-4$.
5. **The classes.** The stored labelling is a complete partition of all 98280 lines, and
   **no class contains a pair at inner product $\pm16$** — this is the whole content of the
   construction, checked class by class, exhaustively.
6. **The axis layer.** $NN^{\mathsf T} = D^2 I$ over the integers, so $Q = N/D$ is a rational
   orthogonal matrix; and for all $1908 \times 1932$ pairs,
   $\big((Qw_i)\cdot w_j\big)^2 \le 48 D^2$ whenever the product is positive — the exact form
   of $\cos \le \sqrt3/2$. The quantity $48D^2$ is checked to fit in a 64-bit integer, so the
   comparison is exact.
7. **The six kinds of pair**, each as an integer inequality, and the final count.

## Data

```
data/golay_basis.txt         12 generators of the binary Golay code used
data/dim14_directions.txt    the 1932 direction vectors, squared norm 8
data/triples.txt             the 644 zero-sum triples, as index triples
data/class_labels.txt        one class label per Leech line, in the canonical ordering
                             produced by the script itself
data/axis_rotation.txt       the integer matrix N
data/axis_rotation_D.txt     the integer D, so that Q = N/D
data/axis_keep.txt           which of the 1932 rotated directions are retained
                             (all of them, since 2026-08-27)
scripts/verify.py            the checker
scripts/axis_rotation.py     how N and D were found: descent on SO(14), then Cayley
scripts/spotcheck.py         an independent numerical spot-check -- see below
```

`scripts/spotcheck.py` does something deliberately different and dumber than `verify.py`:
it materialises actual 38-dimensional unit vectors for a large random sample of the
configuration, over-sampling the pair types that the case analysis treats separately, and
checks every pairwise inner product numerically. It is a guard against a conceptual error in
the case analysis, not a substitute for it — `verify.py` is what proves the theorem, in exact
integer arithmetic. Run it as `python scripts/spotcheck.py [n_sample]`.

The Leech minimal vectors and their canonical line ordering are **regenerated by the
script** from `golay_basis.txt`; nothing about them is taken on trust. The direction set is
Cohn's dimension-14 configuration, reproduced verbatim.

## History of the lower bound in dimension 38

| when | value | who | how |
|---|---|---|---|
| 1982 | 224 608 | Conway–Sloane, *SPLAG* | cross-section of P₄₈ |
| 1998 | 566 652 | Edel–Rains–Sloane [EdRS98] | the three-level chain (n, 8, 2), evaluated with the constant-weight codes of the day — **still the entry in Cohn's table**, and still what Brouwer's page implies |
| Aug 2026 | 570 236 | this project | the same chain with Echols' 2026 A(38,8,8) ≥ 3025 — a +3 584 improvement, now superseded by the row below; see [`../../superseded/dim38-ers-constant-weight/`](../../superseded/dim38-ers-constant-weight/) |
| this work | **591 612** | | the Leech cap construction at codimension 14 |

Dimension 38 is one of two dimensions (with 39) that Echols' own paper passed over, because
he derived his kissing numbers from a fixed level-1 length n₀ = 32 — which for dimension 38
gives 521 084, *below* the record. Dimensions 38 and 39 are exactly the two where the ERS
chain uses n₀ = n. That observation gave 570 236; the construction in this package is
completely different and beats it by a further **+21 376**.

## Context

| dim | previous best | this work |
|---|---|---|
| 38 | 566652 (Cohn's table) / 570236 (this project, ERS + 2026 codes) | **591612** |

For comparison the other constructions in this dimension give: Edel–Rains–Sloane with the
current constant-weight codes, 570236; the 1982 Conway–Sloane cross-section of $P_{48}$,
224608. The same Leech construction was checked in every dimension from 32 to 40 and wins
only in 38 — in 32–37 there are too few triples, and from 39 upwards the code-based records
already exceed the ceiling $3\times196560 = 589680$ of the cap layer.
