# Shared parity repairs improve every dimension from 61 through 68

Checkpoint: 28 October 2026. These constructions replace the preceding
repository constructions in eight dimensions:

| Dimension | Previous bound | New bound | Gain |
|---|---:|---:|---:|
| 61 | 87,056,170 | **88,260,394** | 1,204,224 |
| 62 | 120,620,044 | **121,824,268** | 1,204,224 |
| 63 | 182,537,988 | **184,286,468** | 1,748,480 |
| 64 | 332,170,112 | **335,496,064** | 3,325,952 |
| 65 | 332,238,338 | **335,564,290** | 3,325,952 |
| 66 | 338,302,340 | **341,994,884** | 3,692,544 |
| 67 | 351,420,806 | **355,090,822** | 3,670,016 |
| 68 | 371,552,152 | **375,222,168** | 3,670,016 |

The baseline is commit `816e291`. The total gain across these eight separate
dimensions is 21,841,408 points. Dimensions 49--60 and 69--72 retain their
existing constructions. No unrestricted optimum or fresh worldwide priority
survey is asserted. The manuscript remains its earlier snapshot.

## Reproduction

From the code package directory:

```bash
python verify_shared.py --full-codes
python verify.py
python ../dim64-65-signed-parity/verify.py
```

`shared_layers.py` rebuilds the binary generator, freshly checks all 2^28
parent words, checks every new support and sign code, re-solves the 61D
phase system, and independently recounts the 63D section. Its deterministic
receipt is `shared_layers.json`. `verify_shared.py` runs this in an isolated
temporary copy and compares the entire receipt, along with all inherited
compression inequalities. `--full-codes` additionally rebuilds and enumerates
all **68,719,476,736** words of the [73,36,16] mother code using the independent
scalar C++ enumerator. Its complete histogram must match the supplied one.
The 64--65D verifier separately replays both of its 2^28-word enumerations,
the full 68,224-point owner construction, and this shared-support proof.

The registered 64--65D driver reads this shared receipt; its own older
211-support certificate is retained and checked as a historical construction.
The code package registers 61--63D and 66--68D. There is one shared source of
the new support proof. Both package manifests and their default verifiers
cover the dependency.

## What changes

The original 64D body contains 2^28 dense signs, **30,828** weight-16
supports, 10,416 affine four-supports and 128 coordinate poles. The raw
amplitudes are respectively 1, 2, 4 and 8; all points have squared norm 64
and pairwise inner products at most 32.

The previous augmentation added 211 weight-16 supports. Every new support
met exactly eight old supports in ten coordinates, and no old support could
be reused for another parity repair. That disjointness rule was sufficient
but unnecessarily restrictive.

An extended Hamming [16,11,4] sign code is the kernel of a five-dimensional
dual equivalent to RM(1,4). This dual has the all-one word and thirty
weight-eight checks. Several requested eight-set checks can coexist if
their span still extends to such a dual. An old support can therefore repair
several new incidences while retaining all 2,048 signs.

The new finite search preserves those compatible check spans and solves
the affine phase equations as it fills supports, deletes a few choices,
and refills. The resulting explicit archives are:

| Archive | New supports | Old supports modified | Ten-overlap repairs | Use |
|---|---:|---:|---:|---|
| `global.npz` | 1,835 | 1,792 | 14,680 | 61--65D |
| `expanded66.npz` | 1,803 | 1,792 | 14,424 | 66D |
| `expanded67_68.npz` | 1,792 | 1,792 | 14,336 | 67--68D |

Each archive contains unsigned integer arrays `new_supports`, `new_duals`,
`new_shifts`, `old_supports`, `old_duals`, and `old_shifts`. A support is a
64-bit mask. Its five checks and phase are 16-bit masks in increasing
coordinate order on that support. These fully specify every sign point.
Unmodified old supports retain freely ordered extended Hamming signs.

## Exact support and sign proof

Every exhibited dual has rank five, starts with the all-one word, and has
weight distribution 0, 8, 16. Its sixteen syndrome columns are distinct and
have first bit one, so no one, two or three columns sum to zero. Its kernel
therefore has dimension eleven and minimum distance at least four. The
verifier also enumerates all 2,048 kernel words and checks distance four.
Translating by the supplied phase preserves this distance.

Every new support has size sixteen and meets every other new support in
at most eight positions. Against every one of the 30,828 old supports,
every overlap is at most eight except for exactly eight overlaps of ten.
For each such pair A,B, the verifier derives the unique shared eight-set
H contained in A intersect B and in both sign-code duals. It checks that
the supplied phases have opposite parity on H. Every pair of sign points
then disagrees in at least one common coordinate, giving

    inner product <= 4*(10-2) = 32.

These conditions are checked independently of the search and its indices.
All support intersections are computed twice: packed unsigned population
counts and physical binary-coordinate dot products. The latter are exact
sums of at most 64 binary entries, far below the precision limit of binary32;
there is no angle tolerance. The families require 167,396,040 old/new checks
and 4,911,934 new/new checks in each representation. The expanded-dimension
checks below are also replayed in both representations.

A negative control flips one phase bit in a shared check. The verifier solves
the actual affine sign equations: matching all ten overlap coordinates is
impossible before the corruption and possible afterwards. The resulting
two concrete sign points have inner product **40 > 32**. Thus a phase error
is detected by a geometric counterexample.

Dense and sparse mixed bounds depend only on amplitudes: 1*2*16 = 32,
2*4*4 = 32, and 2*8 < 32. Old/old support distance is unchanged. No old
points are lost in these bodies.

## 61D and 62D

The global family lies in the dual of the self-orthogonal [64,28,16] code.
An exhaustive 635,376-choice four-zero search on the combined support set
finds **10,261** supports avoiding coordinates **(24,42,53,55)**. This is an
optimum for that finite shortening problem only. The dense-code projection
on those coordinates has rank four. Use the same latitude and hexagonal
maps proved in [the preceding continuation](../significant49_72/README.md).
For 61D shorten on coordinates 53 and 55, and use i = c_24 and j = c_42.
For 62D set c_55 = 0 and use the six hexagon labels on (c_24,c_42,c_53),
with c_24 the low bit. In both cases the head is the other 60 coordinates.

For 61D every sparse sign code must have odd total parity, since every dense
restriction to a retained dual support has even parity. The archived phases
are not simply reused. Fix each all-one-check syndrome to one and solve
all surviving shared-check equations simultaneously over F2. The exact
system has **2,769 equations**, rank **2,381**, and is consistent. The replay
constructs a phase for every involved support, checks its odd total parity,
and rechecks every surviving ten-overlap repair. The ordinary unmodified
supports use any single-bit translate. Hence the dense/sparse bounds
14/32 and 14*sqrt(17)/128 remain valid.

In 62D no total-parity condition is needed, so the supplied phases work.
The count identities are

    61: 2^26   + 10,261*2,048 + 8,555*16 + 122 = 88,260,394,
    62: 6*2^24 + 10,261*2,048 + 9,145*16 + 124 = 121,824,268.

## 63D

The explicit `dual_octads.json` lists all 72 weight-eight words in the
parent dual. Membership and distinctness are checked. MacWilliams' identity
applied to the freshly enumerated parent histogram proves that there are
exactly 72, establishing completeness of this small section search.

For each octad H, take the central section perpendicular to its indicator
vector. The dense projection has rank seven and is the entire even-weight
code on eight positions. Exactly 70*2^21 = 146,800,640 dense signs survive.
On a support S let m = |S intersect H|. Odd m contributes zero signs. For
even m < 8, projection is surjective because the sign-code dual has no
nonzero word below weight eight, and the count is

    binomial(m,m/2)*2^(11-m).

For m = 8 the count is 560 if the octad is not a dual check; if it is, the
count is 1,120 for even phase and zero for odd phase. Each unmodified support
can put its octad first, giving 1,120 in that case. This reordering affects
no repaired pair, since all old supports with ten-overlaps are explicitly
modified. The replay enumerates every actual sign word again to check
these projection formulas on the selected section.

The best section retains **1,879,040** added sparse points and all
35,503,776 original sparse points. An independently labelled affine
four-support layer contributes 102,900 points, and 112 ordinary poles
survive. Thus

    146,800,640 + 35,503,776 + 1,879,040 + 102,900 + 112
        = 184,286,468.

## 64D through 68D

In 64D, replacing the 211-support augmentation by the global family gives

    2^28 + (30,828+1,835)*2,048 + 10,416*16 + 128 = 335,496,064.

For 65D choose the dense code independently as the doubly even radical of
the extended BCH code, as in the existing owner proof. This change affects
no dense/sparse inequality. None of the 68,224 owners belongs to the
weight-16 layer, so all survive the new augmentation. Their cap replacement
adds 68,224+2, giving **335,564,290**. The separate 64--65D default verifier
freshly checks that owner argument before accepting this total.

For 66--68D the existing [73,36,16] support code is shortened on respectively
seven, six or five coordinates of its weight-nine dual word. Put the other
64 coordinates first in exactly the order of the rebuilt 64D code, and put
the remaining two, three or four coordinates last. Its 30,828 zero-tail
supports are exactly the old supports already repaired above. There are
3,200, 9,600 and 19,424 additional supports. Every new support lies in the
first 64 positions. The corresponding expanded archive meets every
additional support in at most eight coordinates, checked directly in both
representations. Consequently the old repairs suffice; no additional signs
are removed. The complete Netto four-support layers and poles remain
compatible by the amplitude bounds. The gains are respectively
1,803*2,048, 1,792*2,048 and 1,792*2,048.

## Discovery scope

The optional `../../search/shared_fill.cpp` is the CPU fill/correct discovery
engine. It preserves compatible RM(1,4) check spans, solves binary phase
equations, and tries small deletions and refills. It requires a scratch pool
of records `(uint64 support, uint32 old_indices[8])` plus the original 30,828
packed supports; these large search intermediates are not proof inputs.
Candidate generation enumerated the dual's 3,848,124 weight-16 words and
retained the 112,896 supports with the required eight ten-overlaps. Filtering
against the expanded layers left 62,720, 18,816 and 5,376 candidates in
66D, 67D and 68D. The search seed was 2026102817; finite 90-second walks and
120-second seeded continuations produced the archived families. Replaying
the explicit construction and its proof does not require rediscovering the
search optimum, and no optimum is claimed for these searches.

PackingStar's published fill/correct principle motivated the continuation.
The reusable parity spaces and exact phase systems are specific to this
construction; no reinforcement-learning or PackingStar-code reproduction
claim is made.
