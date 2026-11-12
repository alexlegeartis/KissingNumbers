# Construction-theory pilots for reclaiming dimensions 26–31

12 November 2026. **This package does not improve a kissing-number record.** It records concrete alternative constructions, finite experiments, and small exact certificates that prevent repeating exhausted searches. The successful 28D and 29D exchanges are handled by the campaign's separate verification packages.

`REPORT.md` gives the mathematics, measured results, and scope of every restriction. The key reusable objects are:

- A 928-support signed-octad witness in 31D, giving 218,948 points and therefore falling below the current record.
- An explicit 2,688-point three-core/five-tail layer, together with a proof that its best sampled whole-layer orientation is too expensive against the pinned incumbent.
- The exact 280-head bound for each fixed triangle phase, or 1,680 heads for the six-phase model.
- Six injective matching witnesses showing that transplanting the recovered 28D exchange's core bodies onto every signed `D_k` tail cannot produce a positive gain in those complete supplied pools.

## Short, read-only certificate replays

Run Python without `-O`. Python 3.10+ and NumPy suffice for all verification commands. The commands below are relative to this directory. Apply an affinity that leaves at least one logical CPU available; the two-core examples are suitable for the discovery host.

```bash
OPENBLAS_NUM_THREADS=1 taskset -c 0-1 python verify_small_certificates.py
OPENBLAS_NUM_THREADS=1 taskset -c 0-1 python verify_3plus5_layer.py
```

These commands require no external coordinate files. They regenerate the Golay octads, check all 928 chosen supports and their intersections, verify the Krawtchouk polynomial identity over exact fractions, reconstruct the entire 2,688-point layer, and check its support/sign geometry. They do not claim to replay the source-dependent incumbent deletion costs.

To check those costs and the complete amplification obstructions, provide the pinned source files:

```bash
OPENBLAS_NUM_THREADS=1 taskset -c 0-1 python verify_3plus5_layer.py \
  --with-source --source-dir /path/to/PackingStar

OPENBLAS_NUM_THREADS=1 taskset -c 0-1 python verify_amplification_matching.py \
  --source-dir /path/to/PackingStar \
  --source 28=../../../verifications/improved/dim28-29-november-sextets/data/coordinates28_222120.npz \
  --source 29=../../../verifications/improved/dim28-29-november-sextets/data/coordinates29_231776.npz
```

The source-directory flag may be repeated. A PackingStar repository root, its `25D-31D/new_bounds_20261111` directory, or a directory containing the NPZ files directly is accepted. For differently named or separately located files, use explicit overrides, for example `--source 28=/path/to/file.npz`; the SHA-256 must still agree.

No verification command rewrites evidence. The amplification replay regenerates every advertised candidate from the 28D delta, verifies the matching is injective, and recomputes each matched forbidden product with direct integer arithmetic. All 16 quadratic source points in 26D are supported by exact radical-sign tests. No optimizer or floating-point geometric tolerance is used by these replays.

### Pinned sources

Exact hashes and populations are in `source_manifest.json`. Large incumbent coordinate files are intentionally not copied into this research directory.

| Dimension | Required filename | Population | Provenance |
|---|---|---:|---|
|26|`D26_201180_coordinates.npz`|201180|PackingStar November 11 release|
|27|`D27_209862_coordinates.npz`|209862|PackingStar November 11 release|
|28|`sextet28_222120.npz`|222120|Regained certificate from this campaign|
|29|`sextet29_231776.npz`|231776|Regained certificate from this campaign|
|30|`D30_247992_coordinates.npz`|247992|PackingStar November 11 release|
|31|`D31_251254_coordinates.npz`|251254|PackingStar November 11 release|

PackingStar source repository: <https://github.com/CDM1619/PackingStar/tree/2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6/25D-31D/new_bounds_20261111>.

For the layer-cost proof, an explicit integer flow of value 2,688 saturates the sum of all fiber benefits. The replay checks capacities and conservation directly; no max-flow implementation is required. This proves that no subset of the 21 **whole 128-sign fibers in the exhibited orientation** has positive gain. It does not exclude partial sign fibers or other orientations.

## Reproduce discovery experiments

Discovery writes only to an explicitly supplied output directory. NumPy and SciPy are required; the two C++ searches also need a C++20 compiler. Source-dependent discovery uses the same pinned source arguments as verification.

Regenerate the complete amplification pools and their matchings:

```bash
OPENBLAS_NUM_THREADS=2 taskset -c 0-1 python amplify_delta28.py \
  --source-dir /path/to/PackingStar \
  --source 28=../../../verifications/improved/dim28-29-november-sextets/data/coordinates28_222120.npz \
  --source 29=../../../verifications/improved/dim28-29-november-sextets/data/coordinates29_231776.npz \
  --output-dir /tmp/reclaim-amplification
```

Regenerate fixed-plane triangle pricing and all three-core/five-tail whole-fiber costs:

```bash
OPENBLAS_NUM_THREADS=2 taskset -c 0-1 python triangle_bundle_pilot.py \
  --source-dir /path/to/PackingStar --output-dir /tmp/reclaim-triangles

OPENBLAS_NUM_THREADS=1 taskset -c 0-1 python price_3plus5_fibers.py \
  --source-dir /path/to/PackingStar --output-dir /tmp/reclaim-fiber-prices
```

The disposable triangle candidate/blocker pool is regenerated by the first command and omitted from this package. Its complete pricing statistics are retained in `triangle_bundle_report.json`.

Reproduce the 21-support layer design:

```bash
python build_3plus5_design.py --output /tmp/reclaim-layer-design.json
```

Prepare and run the two bounded C++ discovery pilots:

```bash
python prepare_discovery_inputs.py --support --output-dir /tmp/reclaim-supports
g++ -O3 -std=c++20 support_pilot.cpp -o /tmp/reclaim-support-pilot
taskset -c 0-1 /tmp/reclaim-support-pilot /tmp/reclaim-supports 90

python prepare_discovery_inputs.py --layer --source-dir /path/to/PackingStar \
  --output-dir /tmp/reclaim-layer
g++ -O3 -std=c++20 optimize_3plus5.cpp -o /tmp/reclaim-layer-pilot
taskset -c 0-1 /tmp/reclaim-layer-pilot /tmp/reclaim-layer 90
```

These searches have fixed random seeds but wall-clock stopping rules, so the number of completed moves can vary by host. The saved witnesses are checked independently and do not depend on reproducing the search's timing.

The support preparation obtains the 1,713-support 32D parent and `[32,17,8]` generator from the existing repository path `verifications/improved/dim32-35-signed-layers/source/public`; it does not duplicate those sources here. The shortening audit can be reproduced separately:

```bash
OPENBLAS_NUM_THREADS=1 taskset -c 0-1 python shorten31_probe.py \
  --output /tmp/reclaim-shortening.json
```

## Proof scope

All no-gain statements concern the explicit candidate families or orientations listed in the reports. A matching witness proves that every chosen candidate subset has at least as many incumbent blockers; it is stronger than a search failing to find a gain, but it is not a bound on other candidates. The six-phase triangle bound does not cover arbitrary planes and rotations. The complete 2,688-point layer is valid, but its existence alone does not imply it can be added profitably to PackingStar's incumbent.
