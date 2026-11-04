# Fresh baseline audit, 4 November 2026

No new construction or numerical improvement was found or claimed in this audit. `snapshot.json` records the **existing** 69 registered values in dimensions 25–96, plus the published anchors 47, 48 and 72. The extension below adds all twelve dimensions 5–16 as equally active targets in the same construction language, preserving the earlier data. These are the thresholds a new construction must beat. README/RESULTS and every existing comparison package are preserved.

## Fresh primary-source checks

All eleven workspace HTTPS requests in `retrieval.json` succeeded on 2026-11-04 around 12:36 UTC. Each saved response has a byte length, SHA-256 and timestamp. The browser tool separately returned 502 for Cohn; the successful configured-proxy workspace request provides the saved evidence. No access policy was bypassed.

* [Cohn](https://cohn.mit.edu/kissing-numbers/): all 48 rows parsed. The 22 entries 25–46 now equal this repository's current counts, each citing reference 8 (Kravatskiy's repository). These are listings of these constructions, not an independent tied construction or a proof audit. Dimensions above 48 are absent. Saved HTML: `cohn.html`; parsed rows: `cohn_table.json`.
* [PackingStar](https://github.com/CDM1619/PackingStar/blob/477df34cf2a09bb9476b735a685b538a82da8c55/README.md): latest API commit remains `477df34cf2a09bb9476b735a685b538a82da8c55`, dated 2026-10-28. The 26–31 values remain 201134, 209760, 221864, 231676, 247967, 251169. Thus the existing repository margins remain 2, 47, 124, 21, 5, 19.
* [Lin 32D](https://github.com/QH-Lin/kissing-number-32/tree/f367a9c5a051c0d6821ce71324b9dfa1ead76c2e) and [33–39D](https://github.com/QH-Lin/kissing-numbers-33-39/tree/0c4e2bd486061b3fd55e128ac3497bfb12b73f07): heads unchanged at `f367a9c…` and `0c4e2bd…`. In particular the current 32D README gives 352320, despite the GitHub description retaining the older 352192. The 33–35,37–39 figures remain 370240,390980,415308,499816,605322,768916. Read the actual artifact, not search-result/descriptive snippets.
* [Takhanov–Yun](https://github.com/k-nic/kissing_d25-31/tree/b3768310ba566da5c8947546266c103fbcdce181) and [K. Li](https://github.com/Felpix-Studios/kissing-number/tree/b823640f4e98134b9d55b27288132ee954a8affa): heads unchanged at the repository's cited commits.
* [Brouwer constant-weight table](https://aeb.win.tue.nl/codes/Andw.html) saved fresh. This page is not a complete table at all high-dimensional parameters; absence of a code is not a nonexistence result.

No selected source supersedes a repository value. This is a bounded current-source audit, **not** a complete world-priority certification, especially above dimension 48. Retrieved commit dates do not exceed the context date. Search-engine crawl dates and PDF ingestion dates were not treated as publication dates.

## Targets and comparison hazards

| Priority | Dimensions | Source-backed reason and next useful action |
|---|---|---|
| 1 | 25,26,30,34 | Existing margins over named competing constructions are 1,2,5,10. Require a genuinely larger verified output than 197584,201136,247972,390990. Preserve joint body/head flexibility; a new source release can erase tiny margins. |
| 2 | 27–31 | PackingStar's current README explicitly describes active work on joint removals, midpoint insertion and auxiliary layers. Reflected unions/sextet changes are already included locally. Extend geometry or candidate generation, not repeat those known deltas. |
| 3 | 64–68 | Generic ERS and weight-16 support inputs remain a baseline risk. The old `ers_exposure.py` linear-code tables and prose predate later constructions; do not present their ratios as fresh comprehensive comparisons. Re-evaluate the complete chain with the actual current sign and support parents before judging an improvement. |
| 4 | 88–96 | DGLS is accounted for. In 96D the relevant standard-chain comparison is **141736018112**, not merely its **137438953472**-word sign layer. The current 142162362560 is already 426344448 above that derived chain. A code-only recovery is not a new kissing-number advance. |
| 5 | 73–91 | Cohn's refreshed tail counts update some derived Gamma72 direct-sum floors. In 90D it is now 6218183958 and in 91D 6218188060 (using 8358 in 18D and 12460 in 19D), still below the registered constructions. Historical table comparisons should stay historical. |

## Algebraic inputs and tools beyond the recovered DGLS example

1. **Galand's general coset construction.** [F. Galand, AAECC 2003, pp.235–243, DOI 10.1007/3-540-44828-4_25](https://link.springer.com/chapter/10.1007/3-540-44828-4_25) explicitly gives a method to enlarge binary images of Z/(2^k) codes by adjoining cosets while retaining minimum distance. The primary publisher abstract also says its generalized Kerdock examples are weaker than the best linear codes, while some lifted quadratic-residue examples tie existing linear parameters. This is evidence for a richer assembly operator, **not evidence for a stronger 25–96 input code**. The present `code_search.py` only certifies repeated local translate masks; nonconstant masks, outer codes and code partitions are still planned. Read the full construction before implementing the general operator. Full chapter is subscription content; only publisher abstract/references accessed here. Galand's thesis URL is already linked in the October Z8 review: this is an unimplemented lead, not a newly discovered citation.
2. **Quaternary dual/module families.** [Hammons–Kumar–Calderbank–Sloane–Solé (1994), author-hosted paper](https://neilsloane.com/doc/Me184.pdf), especially §§IV–VI, supplies Z4-dual constructions of Kerdock/Preparata and Delsarte–Goethals/Goethals plus quaternary Reed–Muller structure. A fixed half-degree free cyclic-factor grammar does not cover duals, torsion choices and arbitrary QRM levels. Add those as explicit construction operators with exact distance proofs. Their familiar lengths/distances do not automatically improve a kissing construction: require binary distance >= ceil(length/4), or separately prove the support-layer/section geometry. At 64D a small Kerdock control cannot challenge a 335496064-point baseline.
3. **Ring QR families, with source limits.** [Pless–Qian (1996), DOI 10.1109/18.532906](https://doi.org/10.1109/18.532906), and Galand's references identify extended QR lifts of ring length 32 and 48, which Gray-map into relevant binary lengths 64 and 96. IEEE's page returned a JavaScript robot check, so this audit did not obtain the full theorem/parameter tables. Treat these as retrieval/implementation leads, with no assigned stronger floor until a primary parameter theorem and a verified generator are available.

The strongest immediately actionable algebraic work is a reusable **distance-certified coset/outer-code assembly operator**, followed by exact symmetrized/partial enumerators for large candidates. The earlier campaign already identifies Duursma–Greferath's enumerator paper `math/0302132`; that is not a fresh lead and raw unsigned disjoint-support counts do not implement its character-weighted recurrences. No confirmed missed published code with better eligible parameters was established in this audit.

## Equally active dimensions 5–16

The same [Cohn source](https://cohn.mit.edu/kissing-numbers/), retrieved **2026-11-04T12:36:52.845533+00:00**, supplies the following bounds. SHA-256 of the saved HTML: `94780f35728e349694b5ea119f74e0057b65d5c2d6eb79c9ac0fef9fb3261b3e`. All twelve lower values agree with `common/published.py::COHN` and the earlier campaign's baseline list. No registered improvement in RESULTS exists in this range.

| Dimension | Current lower bound / research threshold | Current upper bound | Exact kissing number known? |
|---|---:|---:|---|
| 5 | 40 | 43 | No |
| 6 | 72 | 77 | No |
| 7 | 126 | 134 | No |
| 8 | 240 | 240 | Yes |
| 9 | 306 | 363 | No |
| 10 | 510 | 553 | No |
| 11 | 604 | 868 | No |
| 12 | 841 | 1355 | No |
| 13 | 1154 | 2064 | No |
| 14 | 1932 | 3174 | No |
| 15 | 2564 | 4853 | No |
| 16 | 4320 | 7320 | No |

These are exact integer **bounds**, not assertions that tau(d) equals the lower bound, except in 8D. Treat 8D as an equally eligible construction/control input with its proved ceiling; an output exceeding 240 would be a verifier or modeling failure. All other dimensions remain active improvement targets, including odd dimensions. Search should use the same operators regardless of dimension; object size and proof availability determine the backend. The former numerical even-dimension survey is not evidence that odd dimensions have been searched.

## One shared missing operator: certified layer gluing

Add **`geometry.layer_glue`** to the common typed construction language. It consumes body objects `B_i` in R^d, rational scales `a_i`, and exact tail vectors `h_i` in R^r, and returns the union of `(a_i*b, h_i)` in R^(d+r). Each body may come from a group orbit, integer configuration, ring/Gray code, quotient coset union, or a previous glued object. No dimension-specific construction engine is needed.

Its proof contract is concrete. For each body certify squared norm `mu_i` and certify `a_i^2*mu_i + ||h_i||^2 = R` for one common positive norm `R`. Supply independently verifiable maximum body inner-product bounds `M_ij`, using distinct pairs when `i=j`. The complete angle obligations are `a_i*a_j*M_ij + <h_i,h_j> <= R/2`, including within-layer pairs. Scales must be nonnegative for this upper-bound formula; negative scales require sign-aware lower/upper envelopes. Check dimensions, exact cardinalities, normalization and all cross-layer obligations. Pair separation then excludes duplicate output points as well. Empty/singleton bodies need explicit vacuous-pair handling.

This recommendation is an inference from the current capability gap, not a new theorem or implemented capability. Small configurations can discharge `M_ij` with exact explicit products; compressed orbit/code objects can discharge the **same obligation** with stabilizers, coset minima, weight enumerators or existing specialized certificates. Refuse an unsupported compressed obligation instead of materializing huge objects or asserting a proof. A first rational/integer backend can avoid new algebraic-number dependencies while retaining a path to exact algebraic scales later.

The operator directly composes low-dimensional latitudes, joint body/head layers, ring-code dense bodies and sparse support families. It generalizes the existing orthogonal union and fixed signed-support recipe. Selecting compatible layers and choosing scales/heads can become construction-program search over this one interface. It offers a useful cross-dimensional next step; it does not by itself promise a new record.

## Replay and remaining limits

Run `python research/discovery_campaign_20261104/baselines/build_snapshot.py`. It parses the saved Cohn table and current RESULTS, verifies all 22 equalities in 25–46, checks all covered direct-sum/DGLS floors remain below registered totals, checks all twelve 5–16D lower bounds against the repository pins, and writes `snapshot.json`. It is a comparison audit, not a fresh execution of the construction verifiers.

No installs, push, publication, modification of comparison packages or existing archive. All writes are under this directory. Source access limitations are ordinary remote content restrictions, not automatic-approval denials; no action is pending. The remaining research blockers are incomplete global source coverage, unavailable full QR/Galand parameter statements, and the absence of a newly verified construction.
