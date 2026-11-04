# High-dimensional coupled discovery

The registered 43D construction now has **2,919,785 points**, improving the
preceding 2,919,782 by three. The generic `search.body_head_exchange` operator
generates all three profitable exchanges from the 30,240 supplied bodies and 305
rational head directions. It receives neither a delta nor a target count.

The separate 43D task discovered the original mechanism and three exchanges.
The engine independently generates the same final point set, performing the
last two steps in the opposite order. This is seeded discovery from imported
body/head/family data. It is not held-out discovery or a worldwide-priority
certificate. The full mathematical replays retain the baseline's P48n
minimum/design inputs and complete-family population identities.

## Positive 43D result

| Stage | Tiny points before | Removed | Added | Gain |
|---|---:|---:|---:|---:|
|Generic one-for-two|3447|1|2|1|
|Generic two-for-three|3448|2|3|1|
|Repeated two-for-three|3449|2|3|1|
|Next bounded search|3450|0|0|0|

The combined delta removes five original tiny points and adds eight. All 118
complete families and five pure heads remain. Each new family shares its
deletion costs. Screening uses floating point and can miss candidates; every
emitted replacement passes exact local integer gates. The engine output remains
an unverified proposal until independent full construction replay.

Both existing representations pass at 2,919,785: the ambient 48-coordinate body
model and the independently reconstructed 40-dimensional coefficient metric.
The latter checks 5,949,525 tiny/tiny pairs, 407,100 tiny/family bounds, 17,250
tiny/pure bounds and 7,621 complete-family/pure cases. Both retain the inherited
population inputs; large shells are not freshly enumerated.

The current point file has SHA256
`042aaff30943205f8d7ce26a4023fb88bcf7db80862484ad21584a1c864135b6`.
The portable package composes the delta against its 3,408-point predecessor:
51 deletions and 93 insertions yield 3,450 tiny points. The former 3,447-point
certificate, delta, reports and notes are preserved byte for byte under
`continuation43/history_2919782/`. `registry43_installation.json` records the
source branch, engine provenance, hashes and inherited-input scope.

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 python3 research/highdim_characters_20261104/engine43.py
OPENBLAS_NUM_THREADS=1 python3 research/highdim_characters_20261104/engine43.py --stage2
OPENBLAS_NUM_THREADS=1 python3 research/highdim_characters_20261104/iterate43.py --extra-rounds 3
OPENBLAS_NUM_THREADS=1 python3 research/43d_head_exchange_20261104/replay.py --delta research/highdim_characters_20261104/engine43_combined_delta.json
OPENBLAS_NUM_THREADS=1 python3 verifications/improved/dim41-46-p48-heads/verify.py
```

The compressed engine receipts contain the complete executable runner programs,
resource profiles, input hashes and proposals. The integration script is a
guarded one-time installer: it refuses to overwrite a changed canonical file.
It is not needed for reproduction after installation.

## Other actual searches

|Dimension|Expanded search family|Result|
|---|---|---|
|33|59 admissible star supports, both leaf parities,456 antipodal candidate lines, shared root-line deletions|370376; gain0|
|34|41 supports, both parities,656 lines, shared root deletions|390990; gain0|
|35|31 supports, both parities,496 lines, shared root deletions|415334; gain0|
|39|13 supports, both parities,208 lines, current deformed dense layer retained|769180; gain0|
|40|91764 distinct algebraic difference proposals from the1134-point auxiliary shell, at norms36,78,84,90|None fits both complete retained large layers|

The 33D–35D and 39D experiments use `search.character_orbits` with the antipodal
sign group. They allow partial stars and both sign parities, rather than the
old mandatory 16-point even-parity stars. Every root family can be removed;
shared deletions are charged once. The solver reaches numerical optimum in
these finite pools. No numerical dual is an exact upper-bound certificate.
Only their root/new-point components are materialized and exactly replayed;
no new full construction is claimed from these zero-gain results.

The 40D proposals are differences `c_i-c_j` of norm 30 auxiliary vectors at
products 12, -9, -12, -15. Exact integer product thresholds screen all 40,590
norm 36 points and, for survivors, the complete 1,084,590-point core. Binary32
products are used only under an explicit absolute integer-accumulation bound
below 2^24. This rejects the specified fixed-core/fixed-medium family; it does
not rule out shared deletions from those layers or other candidate grammars.

```sh
OPENBLAS_NUM_THREADS=1 python3 research/highdim_characters_20261104/root_star_search.py 33 34 35 39
OPENBLAS_NUM_THREADS=1 python3 research/highdim_characters_20261104/probe40.py --inner -15
```

## Comparisons and remaining scope

The strongest available November 4 table capture agrees with the tested
repository counts: 33D 370,376, 34D 390,990, 35D 415,334, 39D 769,180,
40D 1,126,314 and 43D 2,919,782. Its SHA256 is
`94780f35728e349694b5ea119f74e0057b65d5c2d6eb79c9ac0fef9fb3261b3e`.
The primary Lin repository still has pinned head
`0c4e2bd486061b3fd55e128ac3497bfb12b73f07`; its cited counts are lower.
The Yang et al. and Dorofeev–Sun–Wang primary papers were also consulted.
`source_review.json` records the sources and fresh retrieval attempt.

The fresh Cohn request returned HTTP 503; browser retrieval returned 502. Thus
the comparison uses the byte-verified November 4 capture and available primary
sources, with no independent assertion of latest worldwide priority.
The last 43D search found no further proposal in its screened pool. It does
not prove that the construction, finite body/head family or dimension is
optimal. Dimensions 36–38 and 43's separate original research were left to the
parent task; the 43D engine work was added under its explicit integration request.
