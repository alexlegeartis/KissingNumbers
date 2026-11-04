# Executable coverage matrix, dimensions1–24

The matrix has exactly24 rows. Its source is [Cohn's live table](https://cohn.mit.edu/kissing-numbers/), successfully retrieved at **2026-11-04T12:36:52.845533+00:00** through the configured workspace network path. Saved HTML: `../baselines/cohn.html`; SHA-256 **94780f35728e349694b5ea119f74e0057b65d5c2d6eb79c9ac0fef9fb3261b3e**. Browser502/other earlier fetch failures do not invalidate this later successful timestamped response. No new network retrieval was needed for this matrix.

Targets are existing constructions, not new results. Only dimensions1,2,3,4,8,24 have proved exact kissing numbers in this range. The source citations, source hash/time, recipe hashes, proof receipts, provenance and resource limits are recorded per row in `matrix.json`. Separate repository prototypes receive no integrated-coverage credit. Explicitly audited research recipes that execute through the registered runner are eligible; the research allow-list includes `research/coupled13/program_H2.json` plus the17D/18D and Leech orbit/section programs under `research/orbit_controls_20261104/`. Explicit per-recipe resource profiles govern the larger controls, whose receipts use bounded gzip serialization. Planner configuration files without runner nodes and UNVERIFIED forecasts are excluded. The existing CLI `compare` catalogue is intentionally unchanged and may remain partial (for example, a5D recipe can verify exactly while default `compare` reports no recorded5D baseline). This matrix compares replayed constructions against the captured, byte-verified source; it makes no new network retrieval.

| Dimension | Sourced target | Integrated exact target recipe | Provenance | Verification / gap |
|---|---:|---|---|---|
| 1 | 2 | roots01.json, roots_a1.json | generated, imported | exact proof replayed |
| 2 | 6 | roots02.json, roots_a2.json | generated, imported | exact proof replayed |
| 3 | 12 | roots03.json, roots_d3.json | generated, imported | exact proof replayed |
| 4 | 24 | roots04.json, roots_d4.json | generated, imported | exact proof replayed |
| 5 | 40 | gram5_control.json, roots05.json, roots_d5.json | generated, imported | exact proof replayed |
| 6 | 72 | roots06.json, roots_e6.json, section_e8_chain.json | generated, imported | exact proof replayed |
| 7 | 126 | roots07.json, roots_e7.json, section_e8_chain.json | generated, imported | exact proof replayed |
| 8 | 240 | roots08.json, roots_e8.json | generated, imported | exact proof replayed |
| 9 | 306 | support9_search.json | generated | exact proof replayed |
| 10 | 510 | ganzhinov10.json | generated | exact proof replayed |
| 11 | 604 | none | none | prototype_only: current604-point construction and independently replayable proof |
| 12 | 841 | none | none | partial_recipe_only: current841-point construction and independently replayable proof |
| 13 | 1154 | program_H2.json | imported | exact proof replayed |
| 14 | 1932 | ganzhinov14.json | generated | exact proof replayed |
| 15 | 2564 | none | none | prototype_only: 2564-point laminated construction and independently replayable proof |
| 16 | 4320 | stabilizer16.json | generated | exact proof replayed |
| 17 | 5730 | control17_program.json | imported | exact proof replayed |
| 18 | 8358 | control18_program.json | imported | exact proof replayed |
| 19 | 12460 | none | none | prototype_only: current12460 construction; historical11948 mechanism is insufficient |
| 20 | 19448 | none | none | prototype_only: existing19448 layered constructor adapter with exact cross-layer proofs |
| 21 | 31560 | none | none | prototype_only: current31560 construction; historical29768 recipe is insufficient |
| 22 | 49896 | leech_sections_program.json | generated | exact proof replayed |
| 23 | 93150 | leech_sections_program.json | generated | exact proof replayed |
| 24 | 196560 | leech_sections_program.json | generated | exact proof replayed |

`generated` means the recipe constructs a known object; it does not mean discovery of a new construction. `imported` means a certificate/coordinates are supplied as data and verified; this classification is inherited through reachable recipe ancestry. `mixed` combines supplied-coordinate/Gram and generated source ancestors. A numerical search history is recorded separately and is not retroactively turned into autonomous exact recipe discovery. The10D exact recipe is a known-result control with separate numerical line-system rediscovery. The14D exact recipe replays a supplied E7-frame construction. Where a5D Gram import is available, its exactness certifies that data, while the original randomized search provenance remains separately recorded.

Static `roots_a1` through `roots_e8` recipes are literal integer-coordinate inputs at executor level, so their ancestry is classified as imported. The audited inventory also requires exact equality with the independent `root_controls.root_recipe` formula generator. This verified generation lineage is separate from the native `lattice.root_system` routes (`roots01` through `roots08`), which construct roots during program execution. Both routes certify known constructions. The external author's original patch, source commit and author are preserved under `integration/`.

The16D octacode sign recipe produces256 points, below4320; it remains a subtarget recipe even when another16D constructor reaches the target. A registered stabilizer16 control, when present in this matrix, generates4320 points implicitly from affine-flat and quadratic-sign formulas, with symbolic count/separation proof replay rather than a4320-point array. The4D orthogonal-union demonstration imports8 supplied coordinate points, below24. The8D quotient-coset demonstration generates32 points, below240. The exact E8 section chain has target126/rank7 and72/rank6 while retaining8ambient coordinates; output IDs and ambient dimensions are recorded rather than conflated with span rank. Generic integer/Gram import and generic section capability do not imply that every target construction has been encoded. The older19D11948 and21D29768 mechanism notes are below current targets12460 and31560. The registered research `program_H2` route reaches13D1154 by conditional switching, free-zone joining and a coupled lift in Q(sqrt(3)); its supplied coordinate/cap data retain imported ancestry. Current independent replay checks all665281 pairs and exact rank13. It reproduces a known construction. Known17D5730 is reconstructed from Boolean quadratic/RM formulas and verified through the native weighted integer input. The18D8358 construction now uses the native weighted Q(sqrt(3)) backend. A four-action exact orbit recipe generates24D196560; closure, transitivity and196559 seed comparisons cover every source pair. Exact sections supply23D93150 and22D49896. The same signed-pair-normal search reaches only21D27720 and20D17400, below the captured31560 and19448 targets. Maximality of a supplied configuration neither proves that it is globally optimal nor supplies an integrated construction.

Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 python research/discovery_campaign_20261104/coverage_1_24/build_matrix.py --write
OPENBLAS_NUM_THREADS=1 python research/discovery_campaign_20261104/coverage_1_24/build_matrix.py --check
python research/discovery_campaign_20261104/coverage_1_24/build_matrix.py > /tmp/coverage_rows_1_24.json
OPENBLAS_NUM_THREADS=1 python research/discovery_campaign_20261104/coverage_1_24/test_matrix.py
```

`--write` executes eligible currently registered exact controls under their declared resource profiles, separately re-verifies them, writes local receipts and rebuilds the matrix. It refuses to promote missing/unregistered/numerical recipes. `--check` verifies exactly24 ordered dimensions, unchanged source bytes, target/provenance consistency, recipe/evidence/receipt hashes, registry validity and current independent proof replay. Source fingerprints make implementation changes visible. The focused regression checks reject missing rows, a changed source target, false discovery provenance, cross-dimension attribution, omitted integrated runs and an unaudited recipe inventory. Default prints exactly24 row objects as one JSON array after checking. The larger17D/18D and Leech proofs run with explicit profiles; no installations, network requests, pushes or record-registry writes occur.

One shared construction language is used across all dimensions. Resource limits and available proof backends determine feasibility. A gap here describes missing target reproduction, not a theorem of mathematical impossibility or absence of all generic operators. Known controls and benchmark coverage remain distinct from the objective of a new exact improvement.
