# Competitive DGLS coset search, first wave, 4 November 2026

**No new record and no surviving candidate in this bounded run.** This campaign generated candidate families with point-count forecasts above the full 96D repository comparison, then found concrete positive-distance counterexamples for every retained competitive candidate. These finite searches do not exclude other cosets, subcodes, masks, or outer assemblies. Every new file is in this directory; the reviewed engine, prior archive and original generator are unchanged. No install, commit or push was made.

## Competitive objectives and overlap conditions

Let C be the certified free rank-12 Z8 code of length 24, with `|C|=2^36=68719476736` and Gray distance 24. The Gray map sends each ring symbol to an even-parity 4-bit row. The existing second layer is a repeated odd 4-bit mask U applied to every row, so every cross distance between C and C xor U is at least 24 by row parity.

A third complete disjoint compatible parent-sized layer would give **206158430208** points, exceeding the full recorded **142162362560**, without needing compatibility with the existing extra geometric layers. This is a forecast, not evidence that such a layer exists.

Two actual candidate grammars were searched:

1. **Ring cosets `Gray(C+v)`.** Their internal distance is inherited by homogeneous Gray isometry. Different ring cosets are disjoint. Every such even-row layer is compatible with the existing odd-row layer. Cross distance to Gray(C) equals the minimum homogeneous weight of C+v, which must be at least 24. Canonical ring syndromes give `2^36` possible cosets.
2. **Binary translates `Gray(C) xor M` with coordinate-dependent 4-bit row masks.** Internal distance is inherited by translation. Every tested M has a nonconstant row-parity pattern: neither all-even nor all-odd. This proves disjointness from both existing full layers, so all rejection distances are strictly positive and there is no overlap/deduplication ambiguity. Full cross distance varies both parent words; nearest-point screening alone is insufficient.

The ring-coset search sampled uniform canonical syndromes, lifts of the deepest binary Golay residue cosets, and quadratic functions of the coordinate-label bits. It also performed 24 generations of seeded syndrome ascent. The binary-mask search sampled random masks, quadratic block patterns, and masks whose row-parity pattern is a nonconstant Golay residue word.

## Measured first-wave results

| Search | Actual trials or candidates | Outcome |
|---|---:|---|
| Three initial ring-coset families | 12288 | All rejected with actual weight below 24; best 32-information-set upper estimate was 22 |
| Seeded ring-syndrome ascent | 12080 trials, with repeated elites | All rejected; expanded 128-information-set elite upper estimates at most 18 |
| Genuine nonconstant binary masks | 6142 | 5450 rejected by nearest-point witnesses |
| Remaining binary masks, varying both parent words | 692 | All rejected; 21991 pair trials, at most 191 per mask |
| Index-two subcode coset offsets | 88082 across 48 functionals | Generated a bounded compatibility graph using 141979 pair-difference queries |
| Preliminary three-half-coset candidates | 414 | Every one rejected by a stronger 512-information-set screen |

All sampled minima are **upper estimates** of the true minimum distance. A value at least 24 never becomes a distance certificate. Every final rejection instead supplies a specific valid parent message or pair of messages giving a distance below 24.

The exact binary residue covering-radius computation gives radius 4 over all 4096 syndromes. In particular, offsets of the form `4*b` cannot add a complete ring coset: a binary Golay representative within four coordinates gives homogeneous distance at most 16. This explains why a more flexible assembly is needed; it does not rule out general lifted offsets.

## Actual three-half-coset assembly experiment

A partial union of just C and one incompatible ring coset cannot increase the even-layer cardinality: any bad difference word gives a perfect matching between the two ring cosets by addition, so at most |C| points can remain. This matching argument applies to distinct ring cosets with their ring-translation action, not arbitrary partially overlapping binary translates.

The follow-up therefore searched genuinely different assemblies. For a nonzero binary functional ell on the 12 message parities, let `D=ker(ell)`. This is an index-two submodule with `|D|=2^35`. Three mutually compatible D-cosets, together with the full odd-row parent, would give **171798691840** points, above the repository comparison. The run used the 24 coordinate-parity functionals and 24 seeded additional functionals, with random, even, torsion and switching-type offsets. Both D-cosets inside each C-coset were eligible.

The restricted oracle only tests messages satisfying ell(message)=0. Its information-coordinate correction enforces that condition exactly. It generated 414 preliminary triples; each was then tested against 512 fresh information sets on all three pairwise coset differences. Every triple had an explicit bad D-message. This excludes those 414 triples only. No global bound for the subcode grammar is claimed.

A quantitative next candidate grammar is an index-four submodule: five compatible quarter-parent even cosets plus the odd parent would give **154618822656** points. Unlike another whole-parent translate or a two-coset partial union, that forecast can beat the full comparison. It was not run in the first wave; the index-four follow-up below records its subsequent actual search.

## Evidence and reproduce

```sh
OPENBLAS_NUM_THREADS=1 python research/discovery_campaign_20261104/coset96/search.py
OPENBLAS_NUM_THREADS=1 python research/discovery_campaign_20261104/coset96/refine_pairs.py
OPENBLAS_NUM_THREADS=1 python research/discovery_campaign_20261104/coset96/half_cosets.py
OPENBLAS_NUM_THREADS=1 python research/discovery_campaign_20261104/coset96/refine_triples.py
python research/discovery_campaign_20261104/coset96/replay_witnesses.py
```

- `results.json` and `witnesses.npz`: initial constructions, seed, generator/source hashes, covering profile, ascent history and every ring/nearest-point witness.
- `pair_results.json` and `pair_witnesses.npz`: all 692 stronger binary rejection pairs, mask identities, layers and trial counts.
- `half_cosets.json`: all 48 functionals, quantitative counts and all 414 retained preliminary triples, including their complete offsets.
- `triple_results.json`: actual D-message counterexamples for all 414 triples.
- `independent_witness_replay.json`: **PASS**, recomputing 24368 ring-trial witnesses, all 6142 binary-mask rejections and all 414 triple counterexamples. NumPy only reads compressed arrays; direct dot products, Gray words and distances use Python integers. It imports no search code and rejects with explicit exceptions.

The initial search took about 1.5 seconds, pair refinement 0.6 seconds, half-coset generation 14 seconds, full triple refinement 3.8 seconds and independent witness replay 0.9 seconds. The larger elapsed session time includes mathematical setup, implementation and review of the search scope. No giant parent point array is materialized.

There is no competitive candidate to send to a full lower-ring/residue distance certificate in this wave. Remaining work is research, not an environment or authorization blocker.

## Index-four follow-up and reusable generator

The next actual wave tested 118860 offsets across 32 rank-two kernels and 151902 pair differences in retained graphs. A larger run used all eligible offsets for 8 kernels: 29756 offsets and 1006624 pair differences. Its 41 preliminary five-cliques were all removed by actual bad-word witnesses. No frontier survived. The 160837-trial non-torsion ascent then optimized 6 kernels, but fresh 128-information-set queries rejected every collected candidate; the strongest final upper estimates were at most 20. These are finite candidate searches, not exclusions of every index-four assembly.

Pure 4-torsion offsets cannot alone supply 5 quarter cosets: same-parity-label differences have homogeneous distance at most 16 by the exact binary Golay covering radius 4, so the pure family is 4-colorable. Mixed assemblies can still use torsion cosets; the reusable engine retains those vertices and applies only justified same-color edge cuts.

The actual search is now reusable as registered `search.coset_assembly` in `discovery/engines/coset_search.py`, with typed `CosetAssemblyRequest` and explicit `OUTPUT_KINDS=search_result`. It produces compressed unverified candidates and keeps conditional forecasts outside certified facts. See `discovery/docs/COSET_SEARCH.md` and `discovery/recipes/coset96_index4_search.json`. The corrected mixed-family recipe run retains all torsion prototypes, queries 2043 offsets and 8128 pairs, and has no frontier. `index4_full_results.json`, `index4_ascent_results.json` and `index4_witness_check.json` retain the experiment and exact negative evidence.
