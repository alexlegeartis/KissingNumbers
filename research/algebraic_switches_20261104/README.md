# Algebraic discovery in 27D and 31D

The campaign found **209,825 points in 27D** (+18 over the repository) and
**251,196 points in 31D** (+8). These are new explicit constructions relative
to the strongest checked baselines, with complete exact coordinate proofs.
The portable package is
[`dim27-31-algebraic-completion`](../../verifications/improved/dim27-31-algebraic-completion/README.md).

## Positive constructions

The first search transforms whole incumbents with generated Golay sign changes,
signed tail permutations and exact mixed core/tail reflections. It separates
overlap, computes signed cross conflicts and selects individual points. In 31D,
trial 9 exposes four free insertions. The registered
`selection.isometric_exchange` operator independently reproduces them using
the **complete** two-copy union, without a blocker cutoff: 116,504 new image
points, 3,536,678 conflict edges and a checked matching/vertex-cover proof.
It replays both source and output geometry and proves exact rank 31.

The broader signed-support search uses all source six-core/two-tail supports,
all core signs and all signed tail root positions. It breaks symmetry orbits
into individual points. The core oracle is exact support/sign incidence:
only full octad-shaped retained core points can violate an angle; an overlap
of five or six coordinates violates precisely when all overlapping signs
agree. The other core shapes satisfy a checked l1 bound. Explicit bounded
integer products then enumerate every remaining blocker.

|Dimension|New candidate points|Admitted at two blockers|Result|
|---|---:|---:|---|
|27|163,534|92|14 deletions, 32 insertions; **+18**|
|29|1,152,048|0|No proposal in this pool|
|30|2,375,424|32|No gain in this pool|
|31|3,483,640|136|**+4** after the first four insertions|

The portable `lifts27/engine_program.json` and `lifts31/engine_program.json`
execute the registered `selection.blocker_exchange` operator. Inputs are exact
blocker lists and candidate conflicts, with no proposed exchange or target
count. All incumbent points can be deleted; shared deletion costs are charged
once. Oracle geometry is a supplied assumption at this stage, so the output
is explicitly an unverified proposal until full coordinate replay.

Negative-contact sums independently expose the last four 31D insertions.
The resulting eight additions have a compact norm-32 formula:

```
±2(e0+e5+e6+e10) ±2(e8+e20) ±2(e24−e25),
```

with three independent signs. The orthogonal blocks have squared norms 16,
8 and 8. All eight points clear the preceding 251,188-point construction;
no original point is deleted. A different 25-for-29 exchange attains the same
count, but the registered operator prefers the simpler free insertion when
the gain ties. `merged31` verifies the final sorted cube construction.
The four retained points `±4(e8+e20)` and `±4(e24−e25)` complete these eight
additions to a scaled **A3 root system**. Their membership is checked directly
by the portable independent verifier.

These are seeded construction discoveries. They do not load a target count
into the generator and are distinct from the 1–24 held-out benchmark. Golay
actions, root sums and bipartite matching are established mathematics, not
claimed as new mathematical principles.

## Other completed searches

- `images29` and `images30`: twelve whole-image trials each, with cutoff eight;
  no gain. These include signed tail permutations and mixed reflections; 29D
  also moves the normal octad with checked QR affine Golay automorphisms and
  uses tail triality reflections. Multiple-image one-for-two/two-for-three
  checks also found no gain in their stated pools.
- `short_dual40`: all 30,240 shortest dual vectors are reconstructed in the
  exact 40D metric. The 29,106 new vectors may delete points from **all three**
  retained layers. At cutoff 64, 8,694 candidates remain. An exact injective
  blocker assignment excludes gain within that finite pool, independently of
  new/new compatibility. It is not a global bound or an exclusion of other
  alphabets and larger exchanges.
- Four-blocker six/two searches also return no gain: 2,192 admitted points in 27D, 480 in 29D, 2,114 in 30D and 11,752 in 31D. These are numerical finite-pool outcomes, not exact upper bounds.
- Contact closure from the improved 27D and 31D seeds finds no further gain
  in the recorded 64-anchor frontiers. The final odd-split experiments use
  five core and three tail coordinates with explicit body/cost limits; their
  reports state completeness and outcomes separately.

The complete 27D and 31D final proofs cover **53,562,750,010 unordered point
pairs**. Separate direct-int64 calculations check **8,723,928 changed pairs**.
No lattice minimum or unenumerated population is assumed. Original baseline
files and all previous verification packages are preserved.

## Final bounded extensions

`coupled_stars35.py` and `near_star_supports.cpp` release whole middle-support
families together with roots. Both leaf parities are eligible. The first
complete one-middle-blocker catalogue has 1,174 star families and reproduces
415,334 points. The second catalogue has 9,082 admissible families; a stated
6,000-family subset gives 1,933,517 candidate conflicts. Its 90-second MILP
times out without an improving incumbent. This model remains unresolved;
the numerical result is no exclusion. Dense and head layers remain fixed,
and an offending middle support is conservatively deleted in full.

The A3 interpretation also supplies a new generator in `sparse_lifts.py`:
`--core-pair` and `--tail-pair` select four retained orthogonal roots, and the
third root axis runs over every four-coordinate support outside the core
pair. This includes supports absent from the earlier seed catalogue. The
first 31D case tests all 7,315 supports, 234,080 signed core bodies before
pricing, and 75,280 new points; its four-blocker pool has 176 points and no
gain. `a3_batch31.py` explores a bounded nearby set of pair orientations;
its report states the case and wall limits and each actual outcome.

```sh
OPENBLAS_NUM_THREADS=1 python3 research/algebraic_switches_20261104/coupled_stars35.py --middle-limit 2 --out /tmp/new-coupled35
OPENBLAS_NUM_THREADS=1 python3 research/algebraic_switches_20261104/a3_batch31.py --out /tmp/new-a3-batch31
```

The engine now evaluates bounded small integer sums exactly in binary32,
using a proved absolute accumulation bound, while retaining int64/Python
fallbacks. Full-rank elimination stops arithmetic after a full-width minor
but still validates every later input row. Common points of two independently
valid copies skip redundant cross products, with explicit pair accounting.
These changes make the registered full 31D exchange practical without adding
a separate dimension-specific checker.

## Reproduce

```sh
OPENBLAS_NUM_THREADS=1 python3 research/algebraic_switches_20261104/search_images.py 31 --trials 12 --seconds 500 --out /tmp/new-images31
OPENBLAS_NUM_THREADS=1 python3 research/algebraic_switches_20261104/sparse_lifts.py 27 --out /tmp/new-lifts27
OPENBLAS_NUM_THREADS=1 python3 research/algebraic_switches_20261104/run_blocker_pool.py /tmp/new-lifts27
OPENBLAS_NUM_THREADS=1 python3 research/algebraic_switches_20261104/run_engine31.py
OPENBLAS_NUM_THREADS=1 python3 research/algebraic_switches_20261104/search40.py --out /tmp/new-short40
OPENBLAS_NUM_THREADS=1 python3 verifications/improved/dim27-31-algebraic-completion/verify.py
```

`contact_frontier.py` accepts explicit `--source`, `--previous` and `--out`
arguments. `run_engine31.py` uses the stored discovery map, and saves an exact
registered-engine artifact. All source and pool hashes are recorded; stale
PASS receipts are never authority. Research commands write new experiment
directories; full package replay leaves committed evidence unchanged.

The fresh GitHub app lookup confirms PackingStar main remains at
`477df34cf2a09bb9476b735a685b538a82da8c55`, with 209,760/251,169 in 27D/31D.
The pinned November 4 Cohn capture lists 209,807/251,188 and has SHA256
`94780f35728e349694b5ea119f74e0057b65d5c2d6eb79c9ac0fef9fb3261b3e`.
Fresh Cohn retrieval failed. This is a bounded comparison against checked
sources, not an independent latest-worldwide-priority certificate.

The final nearby A3 batch completes all 16 planned cases with no further gain.
It generates 1,787,644 candidate incidences across overlapping catalogues;
this sum is not a distinct-point count. Every admitted pool has complete
explicit blockers, and each numerical selection returns gain zero. These
finite results do not exclude other root pairs, alphabets or larger exchanges.
