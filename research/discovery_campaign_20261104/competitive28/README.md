# Competitive 28D coupled-exchange experiment

No improved construction was found. The incumbent remains **221988** points.
This directory records a performance-directed experiment starting from that
competitive configuration, rather than treating a small Weil orbit count as
success or trying only larger caps.

The current package `verifications/improved/dim27-31-joint-refills` retains a
large Leech core and auxiliary points, and allows coordinated deletions and
insertions through exact cross-conflict graphs. Its supplied sextet script
enumerates only 35 sextets compatible with one fixed normal octad. Repository
Python/Markdown searches found no other full-1771-sextet enumerator; this is a
local history check, not a worldwide novelty claim.

`coupled_search.py` regenerates the Golay octads and all **1771** sextet
partitions from every tetrad, identifying **1736** outside that restricted
35-partition family. It supports direct ordered products of corrected sextet
isometries. Integer numerator matrices are checked against `M*M^T=d^2*I`.
If a composite image is fractional in the source coordinates, both copies are
rescaled to a common exact integer lattice; no fractional image is rounded or
silently skipped. The small regression test explicitly exercises this case.

The actual incumbent's numerator coordinates are all even. Every one of its
35 original sextet images is integral; the fractional-skip rule therefore does
not remove any of those initial 35 branches. The first outside frames and the
four completed experiment words also produced integral images. This is an
observed restriction, not an assertion about every composition or future seed.

For each word, common points are separated exactly. The two exclusive copies
are individually valid by the freshly replayed source proof and exact isometry.
Signed cross conflicts are evaluated using bounded integer products in binary32:
every product and partial dot sum lies below the exact-integer range. Candidates
whose number of blockers exceeds the explicit cutoff are discarded. Surviving
candidates are scanned against every exclusive source point. Maximum matching
and an independently checked equal-size vertex cover then certify the best
compatible subset of this **restricted two-copy union**. An unfinished cross
scan is reported as incomplete, never as a zero-gain obstruction.

`outside_first4_results.json` records two outside frames and both ordered
products. The first two frames commute, so those two products give the same
matrix: four attempts here represent three distinct matrix words. The updated
driver removes identical matrix words before spending further search budget.
All four original conflict pools completed with blocker cutoff 8. Their
restricted candidate counts were 19372, 20170, 19722 and 19722. In every case
the matching covered all candidates, so the finite restricted gain was zero.
Each frame used approximately 1.66–1.70 billion exact cross products and took
12–14 seconds. The complete source geometry was reverified during this run.
`outside_work_cutoff_results.json` records the earlier 16-word run with an
800-million-product per-word cutoff; all those pools were incomplete and
establish no zero-gain conclusion.

Any positive output would require a direct integer check of every inserted
point against every final point, norm/distinctness checks, and source membership
accounting. Unchanged pairs inherit the freshly checked source proof. The script
contains this independent delta checker; the positive six-point regression
exercises it. No positive large output was produced here.

The small phase-cap fix and its reproduction are separate evidence:
`coordinate_orbit_repro.json` shows 13/14 coordinate seeds still unclosed at
cap 4000, while `phase_cap_comparison.json` shows that numerical duplicate block
admission caused avoidable budget waste. Removing that waste does not improve
the 6720-point restricted Weil result or approach the 221988-point incumbent.
`larger_orbit_and_witness.json` further closes coordinate seed 1 at 14560 lines
under cap 20000, but finds a generated vector with product approximately 2/3
against the seed within its first 500 raw images. The whole vector orbit is
therefore numerically inadmissible even if fully materialized. This is a
numerical witness, not an exact algebraic obstruction or a claim about every
other capped seed.

```bash
OPENBLAS_NUM_THREADS=1 python -m unittest discover \
  -s research/discovery_campaign_20261104/competitive28 -p 'test_*.py' -v
OPENBLAS_NUM_THREADS=1 python research/discovery_campaign_20261104/competitive28/coupled_search.py \
  --out /tmp/new-coupled28-run --frames 2 --composites --seconds 120 \
  --blocker-limit 8 --max-checks 4000000000
```

The next structural extensions are to avoid already certified Leech-core/core
cross work through an independently regenerated shell-membership proof, then
couple several image pools or new auxiliary-shell candidates. These steps were
subsequently tested as follows; no generic large-configuration operator is
registered by this experiment.

## Completed multi-image continuation and independent pool obstruction

`multi_image_results.json` records eight distinct matrix words: outside-frame
indices 0, 1, 41, 42 and their cyclic adjacent products. This includes genuinely
noncommuting frames. All eight restricted pools completed. Fresh complete
geometry checks of the incumbent and the generated full 196560-point shell
preceded the search. Exact shell membership allowed **3245528332** certified
core/core cross pairs to be skipped; **9847594496** other exact products were
computed. Discovery took about 111 seconds, excluding those preliminary proof
replays. Every single-image restricted matching had zero gain.

Point-level pooling removed repeated coordinates across the images. It produced
**24178** distinct new candidates, of which 23898 have four blockers and only
56 have two blockers; none has one. The requested two-addition/one-deletion and
three-addition/two-deletion searches completed without gain. Their low-blocker
pool has 56 unique points in 52 buckets; only four candidate-pair products need
checking. This is evidence about that finite pool, not a general small-exchange
ceiling.

The stronger result is an **injective forbidden-owner certificate for all 24178
candidates**, not just the small buckets. The union conflict graph has 97094
edges and admits a matching covering every candidate. The saved witness stores
one distinct incompatible incumbent point per candidate. For any subset of k
added candidates, all k distinct assigned incumbent points must be deleted, so
the gain is at most zero. Candidate/candidate conflicts are irrelevant to this
upper bound and can only make a proposed exchange weaker.

`verify_pool_obstruction.py` independently checks this certificate using direct
integer coordinate products. It trusts no matching algorithm, graph edge table,
or stored PASS receipt. It checks the pinned baseline file, source/candidate
norms, candidate distinctness and disjointness, owner injection, and every
assigned forbidden product. The smallest witness product is **80 > 64**.
`global_blocker_witness.npz` and `global_blocker_matching.json` contain the
witness and discovery summary. Its scope is every subset of these explicit
24178 candidates with arbitrary incumbent deletions; it says nothing about
ungenerated candidates or the points excluded by the blocker cutoff.

```bash
python research/discovery_campaign_20261104/competitive28/verify_pool_obstruction.py
```

This sharpens the structural diagnosis. Of the pooled candidates, 24114 have
zero tail and only 64 change an auxiliary point. Sextet maps fix the tail
coordinates, so they primarily try to reverse existing core deletions. New
auxiliary shells, heights, or exact maps mixing core and tail coordinates are
needed to change that grammar. More optimization on this certified obstructed
pool cannot improve the incumbent.
