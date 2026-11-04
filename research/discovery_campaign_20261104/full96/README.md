# Full 96D replay, not a new record

Both existing independent compressed checkers passed, agreeing on dimensions 88–96 and the full 96D total 142162362560. Layers: dense 137438953472; transversal 4294967296; extra supports 426344448; six-support 2097152; poles 192.

The independent checker took 17.169 seconds and 42632 KiB child peak RSS in the initial run. `resources.json` records that measurement. The subsequent reproducible bridge run records both timings in `bridge_summary.json`.

Run `python -m discovery.benchmarks.run_full96 --outdir /tmp/fresh-full96` from the repository root. It executes both checkers afresh, requires fresh outputs, checks all five layers and complete agreement, and writes source hashes. Comparing saved JSON alone is not mathematical verification. The legacy checkers require Python assertions; optimized mode is refused.

This supplies a full recorded-construction benchmark, not a general mixed-layer operator or a discovery claim. The full 32D 352800 benchmark remains `discovery/recipes/known32.json`.
