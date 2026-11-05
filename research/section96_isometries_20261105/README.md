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
