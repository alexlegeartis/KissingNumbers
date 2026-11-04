# Candidate generation with parity-kernel ring cosets

`search.coset_assembly` runs actual bounded candidate generation. It accepts a compressed free ring generator and binary parity checks on its message coefficients; it generates offsets, deduplicates quotient classes, screens distances, builds compatibility upper graphs and refines proposed cliques. It returns a `search_result` with **empty certified facts** and `Verification('unverified', ...)`. It does not create a spherical code or promote sampled minima to distance guarantees.

The same implementation is available to callers through `search_gray_coset_assemblies(CosetAssemblyRequest(...), context)`. Supported moduli are 4, 8 and 16. Generator rows must have full binary residue rank. Parity checks must be independent. For an r-row parity matrix L and a free rank-k generator G, the search uses the index-`2^r` submodule

```
D = {message * G : L * message = 0 modulo 2}.
```

Its cardinality is `q^k / 2^r`. The requested number of distinct cosets gives a cardinality forecast only. An optional companion cardinality is explicitly an assumption. Parent distance, all coset cross distances, disjointness, and companion compatibility remain proof obligations.

The request provides the generator, parity checks, modulus, target distance, desired coset count and seed. Resource parameters bound generated offsets per family, information-set queries, refinement queries, retained graph vertices, clique attempts and local parity corrections. `RunContext` separately limits dimension, codeword/array allocations, pair work, integer bit counts, trials and elapsed search time. Clique recursion propagates an explicit cutoff. Large correction spaces are sampled without allocating the full population.

Uniform, even, torsion and quadratic coordinate-offset families are available. Each parent ring coset is expanded into all message-parity quotient labels, then deduplicated by its canonical ring syndrome and parity label. For each information set, the oracle matches unconstrained coordinates and corrects an independent set of information coefficients to satisfy every parity equation. Each resulting word belongs to D exactly. A sampled minimum is an **upper estimate** of the true minimum; a witness below the target removes an edge. A proxy above the target produces only an unverified candidate.

The binary residue covering radius is computed exactly when its syndrome space fits the budget. This supplies a narrow torsion-family guard. For offsets whose canonical syndrome is a multiple of `q/2`, two cosets with the same parity label have distance at most `(q/2) * binary_covering_radius`, because `(q/2)` times the residue code lies inside D. Thus an all-torsion assembly has at most `2^r` compatible cosets when that bound is below the target. The implementation can skip an insufficient **torsion-only request**. In mixed requests it retains torsion vertices and applies only justified same-color edge cuts. The family bound never removes torsion vertices from mixed assemblies.

The production handler stores the compressed request, conditional forecast, candidate offsets, sampled upper estimates, exact positive-distance rejection witnesses, statistics and torsion scope. Its verifier checks the known schema and rejection arithmetic but always returns unverified with empty facts. A forged saved exact status does not promote the result; geometric facts placed in `artifact.facts` are rejected. CLI verification exits 2 for these results.

The supplied `discovery/recipes/coset96_index4_search.json` runs the actual five-quarter-coset search around the published DGLS parent. Five quarter-parent even cosets plus a full odd-row companion would give 154618822656 points if certified. This forecast exceeds the repository comparison but is not a record claim. The current bounded recipe run has no surviving five-coset candidate.

Run the recipe normally:

```sh
python -m discovery run discovery/recipes/coset96_index4_search.json --out new_search_receipt.json
python -m discovery verify new_search_receipt.json
```

The second command returns exit code 2 because exact distance remains unproved. `OUTPUT_KINDS` declares `search_result`, allowing the typed planner to generate and retain these search objects without treating forecasts as certified geometry.

`python -m unittest discovery.tests.test_coset_search -v` covers actual generation of a small five-coset assembly (independently checking all 80 resulting points), mixed-request torsion retention, bounded correction sampling, explicit clique cutoffs, malformed/rank/budget requests, conditional companion counts, registration, CLI unverified behavior and forged-proof rejection. The small fixture verifies the implementation; it is not a research claim.
