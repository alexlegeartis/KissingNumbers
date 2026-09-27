# `contributed/` — records sent in from outside

Configurations contributed by people other than the repository's author, kept apart from the four
tiers of [`verifications/`](../verifications/) so that nothing there, and nothing `audit.py` counts,
changes unless the author decides to adopt a record. `RESULTS.md` is not edited by a contribution.

One directory per record, each self-contained: `data/` in the format of the package whose
construction it uses, a verifier that runs from inside the directory, and a README with the count,
how it was found, what was checked and the provenance.

| record | dimension | bound | construction |
|---|---|---|---|
| [`dim27-joint-optimised-two-layer`](dim27-joint-optimised-two-layer/) | 27 | 201 566 | [`dim26-27-iota-triangles`](../verifications/improved/dim26-27-iota-triangles/), re-optimised jointly |
