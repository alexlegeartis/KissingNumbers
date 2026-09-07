# SUPERSEDED — dimension 27 by the two-triangle configuration

**τ(27) ≥ 200 540**, against the published 200 044 — a +496 improvement, held from 2026-08
until 2026-09-08.

**Superseded by this project's own
[`../../improved/dim26-27-iota-triangles/`](../../improved/dim26-27-iota-triangles/), which gives
201 010** — the same template (Leech equator, cap heads on the twelve cuboctahedral directions,
a rotated cuboctahedron on the axis) with the heads rebuilt: the coset classes of two zero-sum
triangles of norm-6 vectors, 2210 heads at squared length 8/3, each removing exactly its
owner, distributed over the four triangles of directions with the side chosen per head.
Better by +470. The partition argument below is still true and is still what fixes the
twelve directions; it is no longer what carries the dimension-27 claim.

---

# A new lower bound for the kissing number in 27 dimensions

This repository is a self-contained verification package for

$$K(27) \ge 200540,$$

improving the previous record of 200044 by 496.

*This directory is one package of the `kissing_verifications` repository and is fully
self-contained; it was written first, as a stand-alone submission, and is left that way.
The construction idea it embodies — that a class of cap heads may be shared by a **zero-sum
triple** of cap directions, so the directions should be partitioned into triples rather than
greedily grouped — is the same one that drives dimensions 73–95; see
[`../dim73-95-gamma72-caps/`](../dim73-95-gamma72-caps/).*

The most important file here is

- `data/dimension27_200540.txt`

and the most important script is

- `scripts/verify_configuration.py`

That script performs an **exact** check of the full 200540-point configuration in
`data/dimension27_200540.txt`. It does **not** rely on the paper's construction, it uses
**no floating-point arithmetic**, and it needs **no third-party packages**. It finishes in
about fifteen seconds. If it passes, the file itself is a valid 200540-point kissing
configuration in $\mathbb{R}^{27}$, and therefore $K(27) \ge 200540$ is proved, even if
some part of the paper were misstated.

## What is in this repository?

```text
.
├── README.md
├── LICENSE
├── CITATION.cff
├── data/
│   ├── SHA256SUMS.txt
│   ├── construction.json
│   └── dimension27_200540.txt
├── paper/
│   ├── kissing27.pdf
│   └── kissing27.tex
└── scripts/
    ├── common.py
    ├── rebuild_from_published.py
    ├── run_all.py
    ├── selftest.py
    ├── verify_configuration.py
    ├── verify_exhaustive.py
    ├── verify_manifest.py
    └── verify_paper_claims.py
```

### `paper/kissing27.pdf`

The write-up, 9 pages. It uses the same normalisation as the scripts.

### `data/dimension27_200540.txt`

The full 200540-point configuration, in the format of Henry Cohn's table of kissing number
bounds: a `Dimension: 27` header, a point count, the inner product coefficients, and then
one comma-separated row of 27 coordinates per point. Rows are **not** unit vectors; each is
normalised by its own length, as in that data set. Every coordinate lies in
$\mathbb{Z}[\sqrt2,\sqrt3,\sqrt6]$ and is written in the same token syntax
(`2*sqrt(3)+sqrt(6)`, `-2-3*sqrt(2)`, …).

### `data/construction.json`

A machine-readable transcription of the small combinatorial data of the paper: the twelve
cap directions and twelve axis points as elements of the ring, the five 496-element classes
of Leech minimal vectors, and the two assignments of classes to directions — the published
one and the improved one.

### `scripts/verify_configuration.py`

The key independent checker. It proves the theorem directly from the data file.

### `scripts/verify_exhaustive.py`

A second, deliberately different proof of the same statement: it assumes nothing at all
about the geometry and compares all $\binom{200540}{2} = 20{,}108{,}045{,}530$ pairs. This
is the only script that needs numpy, and the only slow one (six or seven minutes).

### `scripts/verify_paper_claims.py`

Checks the paper's finite calculations from `construction.json`, and bridges them to the
data file by rebuilding the configuration and comparing.

### `scripts/selftest.py`

Negative controls: fifteen in all. It corrupts the data file in ten specific ways and
confirms `verify_configuration.py` rejects each, then corrupts the data file and
`construction.json` in four more and confirms `verify_paper_claims.py` rejects those.
It also confirms the unmodified file is accepted, which is what makes the other fourteen
mean anything.

### `scripts/rebuild_from_published.py`

Provenance rather than proof: it regenerates `data/dimension27_200540.txt` from Henry
Cohn's published dimension-27 block, showing the file was produced mechanically from that
data and not assembled by hand. It reproduces the shipped file **bit for bit** — the same
SHA-256. It needs numpy and the 97 MB source file, so it is not part of the standard run.

### `scripts/verify_manifest.py`

The SHA-256 of all 17 shipped files of this package, against `data/SHA256SUMS.txt`.
It proves the bytes are the ones the checks above ran on, which is the one thing those
checks cannot establish about themselves. This was a manual step once, and because nothing
ran it, a rebuild of `paper/kissing27.pdf` invalidated the manifest with nothing noticing;
it is part of the standard run now. The manifest covers text files, so an end-of-line
conversion breaks it -- which is what `.gitattributes` is there to prevent.

### `scripts/run_all.py`

Runs the suite.

## Quick start

From the repository root:

```bash
python scripts/run_all.py            # the two checks that need no numpy (~35s)
python scripts/run_all.py --full     # also negative controls and the all-pairs sweep (~11min)
```

or, if you only want the one script that proves the theorem directly from the configuration
file,

```bash
python scripts/verify_configuration.py
```

**Requirements.** Python 3.8 or later (for `math.isqrt`). No third-party packages are needed
for anything except `scripts/verify_exhaustive.py`, which uses numpy, and
`scripts/rebuild_from_published.py`, which uses numpy and the external source file described
below. Do not run the scripts under `python -O`: two of them use assertions as their safety
net and refuse to start if assertions are disabled.

**Checksums.** Every tracked file is listed in `data/SHA256SUMS.txt`. From the repository
root:

```bash
sha256sum -c data/SHA256SUMS.txt
```

The repository sets `* -text` in `.gitattributes`, so no file has its line endings rewritten
on checkout and the hashes hold on every platform.

*Reissued 2026-08-25.* The `README.md` line of the manifest was regenerated on that date,
and only that line. The prose above had described the twelve cap directions as an
icosahedron; they are an A₃ root system, whose cosines are −1, −½, 0, ½. An icosahedron
has no two non-antipodal vertices at −½ and so admits none of the zero-sum triples this
construction is built on — the name was wrong, not the data. `data/construction.json` is
unchanged, as its own manifest line shows, and every number in this package was always
computed from that file rather than from the name.

**The source data set is not included.** `dimensions25-31.txt` is 97 MB and belongs to
H. Cohn's *Table of kissing number bounds*, <https://hdl.handle.net/1721.1/153312>. It is
needed only by `scripts/rebuild_from_published.py`, which takes its path as an argument.
Nothing else here requires it: the package verifies itself from `data/` alone.

## Reading the file

The coordinate syntax is parsed strictly rather than leniently. Every term after the first
must carry an explicit sign, so an unsigned juxtaposition such as `2sqrt(3)` — which one
reader would take for `2 + sqrt(3)` and another for `2*sqrt(3)` — is rejected rather than
silently assigned one of its two meanings; and only ASCII digits are accepted, so a
coordinate spelled with a homoglyph such as U+FF12 FULLWIDTH DIGIT TWO is rejected rather
than read as `2`. Both of these were accepted by an earlier version of the parser. Neither
affects the shipped file, whose 45 tokens are all canonical and are round-tripped through
the parser by `rebuild_from_published.py`; what they affect is whether this package can be
trusted to tell you that some *other* file is the configuration it appears to be. Both are
now negative controls in `selftest.py`.

## Why the exact checker is exact

Every coordinate in the data file lies in $\mathbb{Z}[\sqrt2,\sqrt3,\sqrt6]$; there are 45
distinct coordinate tokens. The scripts represent that ring as integer 4-tuples in the
basis $(1,\sqrt2,\sqrt3,\sqrt6)$, so all ring arithmetic is integer arithmetic. Comparisons
that are not settled by integers alone are settled by certified rational enclosures of the
radicals, computed with `math.isqrt` and tightened until they decide the question.

`verify_configuration.py` first finds the 24 columns that are rational in every row and the
3 that carry radicals, rescales each row to a canonical length, and splits the rows into
three families:

| rows | shape | squared norm |
|---|---|---|
| 194576 | equator $(y, 0)$ with $y\cdot y = 32$ | 32 |
| 5952 | cap $(3y, t)$ with $y\cdot y = 32$, $t\cdot t = 144$ | 432 |
| 12 | axis $(0, a)$ with $a\cdot a = 48$ | 48 |

The rescaling is exact: a row of squared norm $n$ is multiplied by $\sqrt{Nn}/n$ to reach
squared norm $N$, which is rational precisely when $Nn$ is a perfect square, and the result
must be integral. Two rows
that are positive multiples of one another therefore become identical, which is how the
distinctness check catches the same point written at two different scales.

### The heads are Leech minimal vectors

The 196560 distinct vectors $y$ arising as equator heads and cap heads are certified from
the data:

- their shape census is $1104 + 97152 + 98304$, the census of the minimal vectors of
  $\Lambda_{24}$, by shapes $(\pm4^2, 0^{22})$, $(\pm2^8, 0^{16})$ and $(\mp3, \pm1^{23})$;
- the supports of the 97152 vectors of shape $(\pm2^8)$ are **759 octads**, which span a
  $[24,12]$ binary code whose weight distribution is
  $1 + 759x^8 + 2576x^{12} + 759x^{16} + x^{24}$ — that of the extended binary Golay code —
  and which are exactly its weight-8 codewords;
- every one of the 196560 satisfies the congruences that define $\Lambda_{24}$ with respect
  to that code: all coordinates congruent to a common $m$ mod 2, the set of positions with
  $x_i \equiv m+2 \pmod 4$ a codeword, and coordinate sum $\equiv 4m \pmod 8$.

The Golay code is not hard-coded; it is *derived from the file* and then checked to be the
Golay code.

Two classical facts are taken as given, and nothing else:

1. a binary $[24,12]$ code with weight enumerator $1 + 759x^8 + 2576x^{12} + 759x^{16} + x^{24}$
   is permutation-equivalent to the extended binary Golay code, so the congruences above
   define a permuted copy of $\Lambda_{24}$;
2. $\Lambda_{24}$ has minimal norm 4 (Conway and Sloane, *SPLAG*, chapter 10), that is, its
   minimal vectors have squared norm 32 in this scaling.

Given those, the heads lie in a lattice of minimal norm 32, so for two distinct heads $u,v$
the difference $u-v$ is a nonzero lattice vector, $|u-v|^2 \ge 32$, and hence

$$u\cdot v \le 16 .$$

If you would rather not grant even that, run `scripts/verify_exhaustive.py`, which assumes
nothing.

### The six kinds of pair

With that bound in hand, each kind of pair reduces to an integer comparison:

| pair | inner product | condition |
|---|---|---|
| equator–equator | $y\cdot y'/32$ | $y\cdot y' \le 16$ — immediate |
| equator–cap | $y\cdot y'/(16\sqrt6)$ | $y\cdot y' \le \lfloor 8\sqrt6\rfloor = 19$ — immediate given (E) |
| cap–cap | $(9\,y\cdot y' + t\cdot t')/432$ | $9\,y\cdot y' + t\cdot t' \le 216$ |
| equator–axis | $0$ | always |
| cap–axis | $t\cdot a/144$ | $t\cdot a \le 72$ — a $12\times12$ ring table |
| axis–axis | $a\cdot a'/48$ | $a\cdot a' \le 24$ — a $12\times12$ ring table |

Here $\lfloor 8\sqrt6 \rfloor = 19$ because $19^2 = 361 \le 384 < 400 = 20^2$, and condition
(E) — no cap head occurs as an equator point — is checked as a set disjointness.

Only the cap–cap row needs real computation, and only inside the four classes: two caps on a
common direction have $t\cdot t' = 144$, so they need $y\cdot y' \le 8$, that is, cosine at
most $1/4$. That is $4 \times \binom{496}{2} = 491{,}040$ integer dot products of 24-vectors,
a second or two in plain Python. Every other case follows from the $y\cdot y' \le 16$ bound
and the $12\times12$ table of tail dot products.

Between two *different* classes the script also reports the largest $y\cdot y'$, because the
hypothesis of the optimality theorem is that it reaches exactly 16. That search need not be
exhaustive: the Leech certification above already proves 16 is a ceiling, so it stops at the
first pair attaining it. Either way the value returned is the true maximum — the loop stops
early only at a value already proved to be the ceiling.

## What `verify_configuration.py` reports

40 checks, among them:

- the SHA-256 hash of the data file;
- the family split $194576 + 5952 + 12 = 200540$ and the exact norms;
- that all 200540 points are distinct after rescaling;
- the shape census, the 759 octads, the Golay weight distribution and the membership
  congruences;
- condition (E), and that equator and cap heads together are all 196560 minimal vectors;
- 12 cap directions, 496 caps each, 4 distinct classes, pairwise disjoint;
- the maximum head dot product inside a class (8, i.e. cosine exactly $1/4$) and between
  classes (16, i.e. cosine exactly $1/2$);
- $9\,y\cdot y' + t\cdot t' \le 216$ for every pair of directions;
- the two $12\times12$ ring tables, with the largest cap–axis product reported exactly as
  $(24\sqrt3 + 12\sqrt6)/144 = 0.4927\ldots$, against a limit of $1/2$;
- for each of the three thresholds that a configuration could actually violate — cap–cap,
  cap–axis and axis–axis — that it is *not vacuous*, that is, that a value above it is
  arithmetically reachable, so the check can actually fail.

## What the negative controls check

`scripts/selftest.py` takes the real file, confirms it is accepted, and then makes one
targeted corruption at a time — always replacing a row, never adding one, so the file keeps
its declared size and the failure is the defect under test:

1. a repeated point;
2. the same point written at three times the scale (a textual duplicate check would miss it);
3. condition (E) broken — a cap head left on the equator;
4. a cap moved to a direction outside its triple;
5. a foreign vector inserted into a cap class;
6. an equator head of shape $(\pm2^8)$ whose support is not an octad;
7. every axis point moved onto a cap direction;
8. a header whose declared count disagrees with the file;
9. a coordinate written as an unsigned juxtaposition, `2sqrt(3)`;
10. a coordinate written with a non-ASCII digit, U+FF12.

Each must be rejected, and the unmodified file must be accepted.

Control 6 asserts, before it runs, that the support it plants really is not one of the
file's 759 octads — otherwise it would be a control that tests nothing, which is the
failure mode this whole section exists to rule out.

A further four controls run against `verify_paper_claims.py`: `construction.json` with its
axis points forged onto the cap directions, a vector of squared norm 32 that is not in the
Leech lattice placed on the equator, an equator row of squared norm 33, and five duplicated
rows with an inflated header.

### Why these controls exist

Both checkers were put through an adversarial review whose brief was to produce a false
PASS. Both reviews succeeded, and the same mistake caused the worst finding in each: **a
scale factor carried over from a normalisation that no longer applied.**

In `verify_configuration.py`, the cap–axis threshold came from the paper, where axis rows
have squared norm 432, but in this package they have squared norm 48. That gave
$t\cdot a \le 216$ where the truth is $t\cdot a \le 72$. Since $|t\cdot a| \le
\sqrt{144\cdot48} = 48\sqrt3 \approx 83.1 < 216$ always, the check *could never fail*, and
cap–axis went unverified — the tightest condition in the configuration, true maximum
$70.96\ldots$ against a limit of $72$. The review built a file passing all 33 checks while
containing 5952 pairs at inner product $1/\sqrt3 = 0.577\ldots$, by putting each axis point
on a cap direction. That is possible precisely because $\sqrt3$ lies in the coordinate ring,
so $t/\sqrt3$ has squared norm 48 and integral coordinates — and an earlier draft of this
file asserted the opposite, that the ring forbade it.

In `verify_paper_claims.py`, the test for "no axis point is a cap direction" assumed the
factor was 2, from $144 = 4\cdot48$; it is $3\cdot48$, so the factor is $\sqrt3$ and the
loop body could never fire. The same review also showed that the bridge checked the equator
only by cardinality, so a non-Leech vector could be swapped in and still pass; that the
rescaling there had no exactness or integrality guards; and that the row count was never
compared to 200540.

All of it is fixed. The two scripts now share one guarded splitter and one Leech
certification in `common.py`, so they cannot drift apart again; both ring tables assert
that their bound is *attainable* rather than vacuous, which is the check that would have
caught the original error; and every exploit the reviews produced is now a permanent
control above.

## The construction, in one paragraph

The configuration has the shape used by the record configurations in dimensions 25 through
31: an equator of Leech minimal vectors $(x,0)$, caps $(\sqrt{2/3}\,u, \sqrt{1/3}\,\hat w)$
with $\hat w$ a unit vector in $\mathbb{R}^3$, and axis points $(0,\omega)$. Each cap
direction $\hat w$ carries a set $U_{\hat w}$ of 496 Leech vectors, and every vector used as
a cap head must be deleted from the equator. Two directions can share a set only when their
cosine is at most $-1/2$, and at most three unit vectors are pairwise at cosine $\le -1/2$,
so a set serves at most three directions. The published configuration groups its twelve
directions as $\{3,3,2,2,2\}$, spending five sets and $5 \times 496 = 2480$ deletions. But
the twelve directions admit exactly two partitions into four triples, so four sets suffice
and only $4 \times 496 = 1984$ deletions are needed. The 496 vectors freed by dropping the
fifth set return to the equator, and

$$194080 + 496 \;+\; 5952 \;+\; 12 \;=\; 200540 .$$

No new spherical code is constructed: the Leech vectors, the twelve cap directions, the
twelve axis points and four of the five classes are reused verbatim from the published
dimension-27 configuration. Only the assignment of classes to directions changes.

## History of the lower bound in dimension 27

Every configuration in this family satisfies the master formula

$$τ(24+k) \;=\; 196560 \;+\; 2 \cdot (\#\ \mathrm{triples}) \cdot s \;+\; τ(k),$$

with $s$ the class size, so the history is the history of two numbers: how large a class
anyone could find, and how well the cap directions were grouped.

| when | class $s$ | grouping of the 12 directions | τ(27) | who |
|---|---|---|---|---|
| 1971 | — | — | 196560 | Leech–Sloane [LS71], inherited from Λ₂₄ |
| 2011 | 480 | — | ≈ 199 000 | Cohn–Jiao–Kumar–Torquato |
| 2016 | 488 | — | ≈ 199 400 | Kallal–Kan–Wang [KKW16], arXiv:1608.07270 |
| Nov 2025 | **496** | {3,3,2,2,2} — five classes, 7 gain-units | **200044** | Ma et al. [Ma25], arXiv:2511.13391 ("PackingStar") — the current table entry |
| this work | 496, unchanged | **{3,3,3,3} — four classes, 8 gain-units** | **200540** | the twelve directions, an A₃ root system, admit exactly two partitions into four triples |

The 2011 and 2016 rows are approximate because those papers report class sizes rather than
a grouping; what is certain is that 200044 is the published value and that it uses five
classes where four suffice.

**This is the last slack in the family.** By the same counting the gain is at most
$\lfloor 2\,|W|/3 
floor = \lfloor 2\cdot 12/3
floor = 8$ classes' worth, now attained,
so dimension 27 is closed to this argument. Running the same accounting across dimensions
25–31 reproduces 198550, 204520, 209496, 220440 and 238350 *exactly* and finds slack only
in 25 (+2, a different mechanism — see [`../dim25-cap-level/`](../dim25-cap-level/)) and
here.

## Provenance

The configuration is not new coordinates. Running

```bash
python scripts/rebuild_from_published.py path/to/dimensions25-31.txt /tmp/out.txt
```

against the dimension-27 block of Cohn's data set reproduces `data/dimension27_200540.txt`
byte for byte, with the same SHA-256. The Leech minimal vectors, the twelve cap directions,
the twelve axis points and four of the five classes all come from that block unchanged; the
script only re-assigns classes to directions and returns the fifth class to the equator.

## Checksums

The SHA-256 checksums are in `data/SHA256SUMS.txt`. For the main configuration file:

```text
3412478bd6dc8c105378f0047d13da2e757c1bf3a180feac4407613cefc387e2  data/dimension27_200540.txt
```

## Licence and citation

MIT (see `LICENSE`). `CITATION.cff` carries the metadata GitHub uses for its "Cite this
repository" button.

## Bottom line

If you only want the shortest path to confidence in the theorem, run

```bash
python scripts/verify_configuration.py
```

If that script succeeds, then `data/dimension27_200540.txt` is an exact valid kissing
configuration with 200540 points in dimension 27, and

$$K(27) \ge 200540$$

follows immediately.
