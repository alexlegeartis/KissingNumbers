# Sol continuation, 4 November 2026

**No new certified SOTA was found.** The unresolved objective is a stronger
construction, not an access or permission issue.

The original extension was restored into an isolated local worktree from
commit `e25c2924f377a1061b6c60429c950d34df7515e9`. The Library archive SHA-256
is `fb8aec51aea76f92306939aa082f562d0c42d0ba6dac2bc7d7ac0b83605f124a`.
All 405 manifest files passed size/hash checks before restoration. The original
archive and the original KissingNumbers checkout remain preserved outside this
worktree. The restored suite passed 168 tests; the repository audit passed all
70 existing registered claims.

The continuation adds three compact construction/search operators, taking the
registered count to 32, and a separate target-blind 1–24 benchmark. There were
no subagents, software installations, Actions changes or PRs. GitHub publication
was subsequently explicitly authorized by Alex.

## Held-out search and exact validation

`blind1_24/search.json` records 200 actual construction programs and 24 ordered
dimension results, using no comparison inputs or incumbent coordinates. Search
took about 187 seconds of summed execution/replay time. `evaluation.json` was
produced by a separate evaluation phase that replayed every saved artifact and
then read the SHA-verified captured table.

- Combined exact matches: 1–9 and 16. Some are formula controls.
- Actual support-optimizer matches: 4,7,8,9, including 9D306.
- Historical 12D840 is below the current 841 target.
- The remaining rows are visible failures; 24D7040 is far below 196560.
- Integrated target *replay* coverage, after adding the 9D search recipe, is
  1–10,13,14,16. This broader replay number is not the held-out discovery score.

All 179 package tests pass. Eleven new focused cases also pass under Python
`-O`. Source and candidate geometry remain independently replayable. Exact
validity, historical provenance, target attainment and novelty remain separate.

## Search priorities from neighboring bounds

The diagnostic uses all 48 captured lower/upper pairs and exact rational growth
ratios. It handles the exact dimensions' merged table cells correctly. The
one-step gap-expansion order begins 31,37,29,26,28,27,33,34,35,25.

For 31D, lower growth is `251188/247972 = 1.0129692062`, while upper growth is
`2213487/1537707 = 1.4394725393`, giving gap expansion `1.4210427430`.
For 37D the corresponding factors are `505503/501687 = 1.0076063362` and
`18498316/13017098 = 1.4210783387`, giving `1.4103507368`.
The evaluator also saves centered two-step diagnostics so a single plateau
need not be interpreted in isolation. These are research heuristics only.

Source: `research/discovery_campaign_20261104/baselines/cohn.html`, captured
2026-11-04 12:36 UTC, SHA-256
`94780f35728e349694b5ea119f74e0057b65d5c2d6eb79c9ac0fef9fb3261b3e`.
Fresh requests failed with HTTP 502 and 503; `retrieval_status.json` records the
fallback. The diagnostic does not claim a newly fetched global leaderboard.

## Competitive construction experiments

Repository history and the restored campaign reports were consulted first.
The 28D sextet pool already has a complete injection obstruction and fixes
tail coordinates. The 14D two-copy campaign and heterogeneous-weight barriers,
the nonclassical 13D seeds, the 17D layered ceiling and higher-dimensional ring
obstructions were read before choosing this continuation. Text searches found
no existing general contact-sum operator; this local check is not a literature
novelty claim. The 13D layer counts and separate sign-code ceiling were already
known in KNOWLEDGE section 66, so generating another sign code alone would not
solve the coupled construction problem.

| Experiment | Actual generated/admitted pool | Outcome and scope |
|---|---:|---|
| 13D contact completion, nonclassical seed `(0,2)` | 8552 new points; all admitted | Complete joint graph has 200978 candidate conflicts. All original components are removable. A 60-second full-pool MILP returned gain 0 and a loose numerical gain bound 771; unresolved, not an exclusion. |
| 14D registered contact-exchange DAG | 3000 new points, plus 1932 originals | Exact graph has 45130 conflicts. Numerical optimum returned 1932. Generation stopped at its explicit 3000-point cap; no unrestricted conclusion. |
| 31D core–tail reflection/contact pool | 6862 new points, all retained in full run | Exact injective old-owner certificate excludes gain from every subset of this explicit pool, with arbitrary incumbent deletions. No conclusion about other points. |

The 31D pool comprises 4929 points from 24 exact Householder maps mixing three
Leech coordinates with one tail coordinate, and 1933 auxiliary contact sums.
For a four-entry sign normal v, `M=2I-vv^T` satisfies `MM^T=4I`; the existing
even coordinates make the sampled images integral. This changes auxiliary
coordinates, beyond the old tail-fixing sextet grammar. Every generated point
has norm squared 128. A full source replay verifies the existing 251188-point
incumbent and covers 31547580078 unordered source pairs.

`mixed31_full/witness.npz` assigns each of the 6862 candidates a distinct
incompatible incumbent point. `contact_campaign.py --mode replay` checks the
source/witness hashes, integer shapes, common norms, new-point distinctness and
every assigned forbidden product directly. It calls no graph or matching
solver. This proves any k insertions require at least k old deletions.
`mixed31_full_replay.log` records successful replay under optimized Python.
The earlier 178-point cutoff run is retained as a narrower historical step.

The 13D pool has 24 axis points, 4856 other signed four-support points and
3672 full-sign fixed-height points. It uses the freshly verified nonclassical
1154-point seed; neither antipodality nor an immutable equator is imposed.
The earlier 3000-point solve and an aborted pre-budget run are preserved.
The aborted run ended before optimization because its 8552-point pool exceeded
the initial 5000-point guard; its witness is not a negative certificate.

A further full-pool 13D run adds exact Hamming-ball clique cuts to the numerical
model. A fixed-height sign layer has dot `norm-2*HammingDistance`; a radius-one
ball has pair distance at most two and is a conflict clique for norm > 8.
This improves the model without fixing any original point. Its separate
`contact13_cliques/report.json` records 4096 such cuts; the run again returned
gain zero at its 60-second limit, with numerical gain bound 771, so this
strengthening did not resolve the pool.

## Reproduction

The supplied 32D projection observation was independently audited with the
standard-library `audit32_projection.py`, including under Python `-O`.
It reconstructs the dual by binary elimination, enumerates both code and dual,
checks the entire dual weight histogram and the exact support census, and
independently counts blockers over all 131072 parent words and all 4096 sign
patterns on examples of every shortening type. All stated counts pass:
rank12 gives 416 blockers; rank11 weight8/10/12 gives respectively
320/512, 192/640, 64/768; rank10 gives 128/512. Exactly 264 twelve-coordinate
supports have projection rank10. The minimum of64 applies to equal-amplitude
norm32 vectors on twelve coordinates with the complete dense layer retained.

The radius-one claim through support size14 follows from dual minimum8.
Plotkin's weight sum excludes shortening dimension4. At dimension3, equality
forces size14, seven nonzero words of weight8, and each nonzero syndrome column
twice. At dimension2 all three nonzero column types must occur, since an absent
type would force two disjoint weight-at-least8 words on at most14 coordinates.
Dimensions0/1 are immediate. Thus every syndrome is zero or a single column,
so the projected code has covering radius at most1. This gives the general
fixed-layer obstruction `sum(abs(x))-2*max(abs(x))>16` and the equal-amplitude
support12–14 corollary. Neither argument excludes a coupled exchange.

The earlier16-coordinate radius-two obstruction and the dense code's covering
radius6,35-optimality and single-coset injection were already in KNOWLEDGE
sections130 and132. This audit refines fixed-layer information and supplies no
new construction or literature-priority claim. A future generic projection
backend should use the existing syndrome machinery and expose blocker prices;
it must not discard a candidate merely because a mutable dense layer blocks it.

The follow-up affine-kernel claim also passes an independent complete audit.
The installed C++ compiler builds a temporary exact-rank checker; no software
is installed. All 6420736 pairs of the 3584 six-dimensional blocker kernels
give intersection dimensions `{0:3785280,1:1956448,2:608448,3:69360,4:1200}`.
For any two affine blocker fibers the intersection is empty or a translate of
the direction-kernel intersection, hence has at most16 words. Two64-word
blocker sets therefore have union size at least112. Combined with the preceding
complete blocker classification, fewer than112 dense deletions can unlock at
most one point from this equal-amplitude weight12 alphabet. This is a bounded
exchange obstruction, with larger exchanges and other alphabets still open.
`audit32_kernels.json` records the complete pair census and source hashes.
The 1200 pairs with intersection dimension4 also have 1200 distinct canonical
eight-dimensional spans. Hence any eight-dimensional space contains at most
two of these kernels: three would produce repeated pair spans. An affine
eight-flat contains at most four cosets per such six-dimensional kernel,
so deleting its256 dense words can free at most eight minimum64-blocker
weight12 candidates. This says nothing about arbitrary256-word deletions or
other blocker costs; it is not a construction-wide ceiling.
The additional octad-layer cost check passes too. The certificate has1704
complete even-sign octad supports and the same dense code in the same
coordinate order. Direct even-sign enumeration gives minimum blocker
multiplicities4,2,8,1 at overlaps5,6,7,8. Across all3584 dual weight12 supports,
at least140 octads overlap in five or more positions, and the weighted sparse
blocker price is at least498. The disjoint dense and complete-octad layers
therefore require at least562 deletions for one such insertion. This prices
individual points; it does not bound shared-deletion gains for a group.

## Reproduction commands

```sh
OPENBLAS_NUM_THREADS=1 python research/sol_discovery_20261104/contact_campaign.py \
  --mode 13 --out /tmp/new-contact13 --seconds 60 --pool-limit 12000
OPENBLAS_NUM_THREADS=1 python research/sol_discovery_20261104/contact_campaign.py \
  --mode 13 --out /tmp/new-contact13-cuts --seconds 60 --pool-limit 12000 --clique-cuts
OPENBLAS_NUM_THREADS=1 python research/sol_discovery_20261104/contact_campaign.py \
  --mode 31 --out /tmp/new-mixed31 --cutoff 2048
OPENBLAS_NUM_THREADS=1 python -O research/sol_discovery_20261104/contact_campaign.py \
  --mode replay --out /tmp/new-mixed31
python -O research/sol_discovery_20261104/audit32_projection.py
python -O research/sol_discovery_20261104/audit32_kernels.py
OPENBLAS_NUM_THREADS=1 python -m discovery run \
  research/sol_discovery_20261104/contact14_program.json --out /tmp/new-contact14.json \
  --budget max_pairs=30000000 --budget solver_seconds=20
```

The first 31D pass took 213 seconds using exact integer matrix products. The
same pool's complete follow-up took 53 seconds using explicitly bounded exact
binary32 integer products and then direct integer witness replay. The complete
8552-point 13D graph/solve took 67 seconds. These are bounded experiments, not
an exhaustive search across constructions or a proof that existing bounds are
optimal. Further progress needs a construction grammar that can escape the
observed owner obstructions and optimize difficult mixed components.
