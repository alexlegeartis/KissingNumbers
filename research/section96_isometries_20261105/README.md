# Section-first search for 90D and 92D

This is a **conditional search direction**, not a new bound. Let C be the
2^36-word even Gray image of the free Z8 Golay code in the DGLS construction.
For m=4 or 6, retain the words whose first m four-bit blocks have weight 2.
Their sign images lie in the (96-m)-dimensional subspace where those block
coordinate sums vanish. The triangular generator makes the first m ring
coordinates independent and uniform, so this section contains 2^36(3/4)^m
words. Each word has squared norm 96 and internal distance at least 24.

Choose a blockwise isometry T of the even four-bit words, allowing the last
blocks also to map to odd words. Require an even image for each of the first m
blocks. The section of T(C) then has the same count and internal distance.
If and only if *every* pair between the two sections has distance at least 24,
their union is a kissing configuration of twice that size. Thus the conditional
counts are 43,486,543,872 in 92D and 24,461,180,928 in 90D. The respective
repository comparison figures are 21,880,490,288 and 12,306,670,328. Any result
must be checked against current literature before it is called a record.

The crucial relaxation is that cross-distance is required **only between words
in these sections**, not across all of C and T(C). `search.py` shares one block
map among coordinates congruent modulo the requested group count. It proposes
maps by mixed-integer optimization and adds a valid cut from each explicit
forbidden pair. The cut accounts for the possibility that the candidate word
leaves the section after a map changes. `oracle.cpp` uses randomized local
search with first-m coordinates initialized inside the section. It does not
exhaust all pairs or all map choices. A missing witness, solver status, or
numerical candidate is not a certificate; a separate exact cross-distance
proof and full construction check are required.

From this directory, compile and run, for example:

    g++ -O3 -std=c++17 oracle.cpp -o oracle
    OPENBLAS_NUM_THREADS=1 python -u search.py 4 4
    OPENBLAS_NUM_THREADS=1 python -u search.py 6 6

The first parameter is the number of periodic coordinate groups and the second
is m. Unresolved proposals are saved as NumPy arrays solely for further tests.
The source currently uses a fixed number of randomized oracle restarts; no
result is promoted automatically. All files are local and no leaderboard is
modified by running this experiment.

For the 86D section (`m=10`), the same conditional count is 7,739,670,528,
compared with the repository figure 6,225,165,512. Compile the additional
`oracle_m10.cpp` as `oracle_m10` before running `search.py 10 10`. For each
randomly chosen pair of valid ten-coordinate prefixes, this oracle exhausts
all 64×64 combinations of the remaining message coefficients. The prefix
pairs themselves are sampled, so failure to find a forbidden pair does not
certify global compatibility. This extra search targets 86D directly rather
than requiring a valid full 96D union.

An optional third argument `s` restricts the transformed layer to the subcode
where message coefficients m through m+s-1 are even. The original section is
left intact. Because the first m output coordinates do not depend on these
later coefficients, the transformed section has exactly 2^(-s) of the full
section's size. With s=2, the conditional counts are 27,179,089,920 in 92D
and 15,288,238,080 in 90D. A compatible map is still required. Compile
`oracle_index.cpp` to `oracle_index` and run `search.py 4 4 2` or
`search.py 6 6 2`; the randomized oracle remains incomplete. In particular,
absence of a witness is not a lower-bound proof.

A fourth optional integer argument is a bit mask selecting which `s` message
coefficients are even, so long as all selected indices are at least m. For
example, `search.py 4 4 6 4032` selects coefficients 6–11. Different masks
specify different subcodes; all have the same conditional section cardinality
under the triangular-generator argument. For s=6 the corresponding potential
counts are 22,083,010,560 in 92D and 12,421,693,440 in 90D, still contingent
on finding and exactly verifying cross-compatibility.

## A common symmetry for partial-layer exchanges

There is a common subgroup even when T differs from block to block. The 4C
subcode has 2^12 words. Adding a word of 4C complements selected four-bit
blocks, and every allowed block isometry commutes with complementing. It also
preserves each weight-two section. Consequently 4C acts freely on both the
original section and every transformed even-coefficient subcode section.

In 92D the original section has 5,308,416 orbits of 4,096 points; an index-64
transformed section has 82,944 such orbits. In 90D the respective numbers are
2,985,984 and 46,656. This suggests a search that retains or replaces *parts*
of both layers in complete common orbits rather than requiring the entire
transformed subcode to coexist with the entire original section. These counts
and equivariance are verified by `python verify_common_orbits.py`. They do not
establish which orbit pairs are compatible, a successful exchange, or a bound.
