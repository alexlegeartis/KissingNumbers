# Independent exhaustive 30D audit

Run from any working directory:

```sh
python verifications/improved/dim29-30-native-partial-design/independent30/audit.py
```

`audit.py` loads only the package's shared `../data/coordinates30_248180.npz`, `source30_248084.npz`, and `delta30_248180.npz`. It divides the even candidate and source coordinates by 2, to squared norm 32; delta coordinates already have this scale. No large arrays are duplicated here.

It verifies cardinality, uniqueness, every norm, the exact source exchange (384 added, 288 removed), and all 30,796,532,110 unordered candidate pairs. The maximum distinct-pair inner product must be at most 16. Full pair multiplication uses binary64 as exact integer arithmetic: input bounds are checked; squared norms are 32, so Cauchy–Schwarz bounds each sum of absolute products and every partial sum by 32. There is no numerical tolerance. An independent int64 replay checks all added-versus-retained and added-versus-added pairs.

`audit_result.json` and `audit.log` record the latest completed execution, including SHA-256 of all three input files. The script writes its receipt alongside itself. A PASS establishes feasibility and the exact +96 exchange; it does not by itself establish literature priority.
