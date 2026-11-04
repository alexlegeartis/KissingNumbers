# Contact-guided 12D search

**No 842-point kissing arrangement was found. The verified lower bound remains
841.** This campaign supplies a slightly better numerical near-miss, an exact
local exclusion certificate, and completed discovery/recovery controls.

## A small numerical improvement, with an exact bound

The supplied published 842 realization has maximum cosine
`0.5009010186022465`. Smooth minimax continuation reduced this to approximately
`0.500900474851737`. Every pair in `evidence/improved842.npy` is independently
certified, using integers on normalized dyadic rays, to have cosine at most

```
500900475 / 1000000000 = 0.500900475.
```

This closes only about 0.06 percent of the gap to 0.5. The same exact check
finds **10,452 pairs still above 0.5**. It is not a new kissing lower bound or a
claim to improve every independently available 842 near-miss.
See `improved842_certificate.json`.

## A rigorous local obstruction

The stress fit was converted into rational pair weights and checked with
integer interval arithmetic. Around the literal stored rays of the published
near-miss, any valid 842-point configuration must move at least one matched
point by more than **0.024484 in chord distance**, about 1.4 degrees.

The proof does not assume the numerical stress is exactly balanced. It includes
its entire residual and a rigorous quadratic remainder. It excludes a specified
neighborhood; it does not prove global impossibility, global optimality, or a
necessary contact-graph change. The corresponding certificate around the
improved near-miss has radius 0.024409.

[LOCAL_EXCLUSION.md](LOCAL_EXCLUSION.md) gives the derivation and independent
replay commands. Both small input files are kept in `evidence/`, with their
hashes and exact rational bounds in the two `*local_certificate.json` files.

## Discovery and recovery controls

The original PyTorch search functions were run in float64 with their complete
35,000-step default search schedule. For the paper's batch-64 scaling, the
first 9,000 retained steps ran on the RTX A4000. The saved raw coordinates,
Adam moments, and step counter were then resumed as 64 independent CPU
trajectories, each dividing its loss gradient by 64. A separate check reproduced
batched-versus-split Adam updates exactly on a small CPU example. Switching
GPU/CPU changes roundoff and is not a bitwise GPU replay.

Those 64 starts finished with a best maximum cosine about **0.52592131**.
Another eight trajectories, using the released script's default batch-512
gradient scaling, completed the same schedule on the GPU and reached about
**0.53371943**. Batching only couples these losses through that scalar divisor,
so the eight trajectories are a subset-sized test of the 512 scaling, not 512
starts. Neither run produced an output below the released save threshold
`0.501`; consequently there was no eligible candidate for its separate
750,000-step polishing script.

This is a failed reproduction in the tested starts, not a disproof of the paper
or evidence that its much larger random search could never succeed. The paper's
841 coordinates remain independently certified. These controls do not establish
a reliable discovery success rate for our environment.

Separate tests used direct feasibility optimization (128 starts) and the prior
projected-Adam explorer (64 starts) for 841, without a discovery success.
Recovery from the known 841 was more successful: **5 of 16 perturbed starts**
returned configurations that pass the exact kissing test. One successful start
used coordinate noise with standard deviation 0.03. This verifies local recovery
capability; it does not demonstrate discovery from an 840 construction.
An exactly checked recovery is stored as `evidence/recovered841.npy`.

## Changes to contacts and bridges

The full minimax linear programs did not return solutions within their requested
budgets. We instead fitted nonnegative weights on near-maximal pairs by sparse
least squares. Their residuals and complementarity gaps are recorded explicitly.
Three different contact bands consistently identified the same most stressed
points. The exact exclusion proof is valid independently of whether this fit
is optimal.

The continuous searches used these weights to select regions and contacts:

| Search | 842 starts | Best maximum cosine |
|---|---:|---:|
| Rebuild 16–128 stressed points using soft motions, rotations, noise, or splitting | 192 | about 0.50094157 |
| Insert a point into the known 841, release all blockers and a growing boundary | 192 | about 0.50125320 |
| Force 1–8 stressed contacts open, then release the added constraints | 144 | about 0.50094111 |
| Favor fewer blockers, starting from five exactly verified recovery variants | 96 | about 0.50119611 |
| Release low-conflict algebraic trade candidates continuously | 64 | about 0.50548060 |

There were **688 new 842 construction trials** in these batches. Most moderate
patch changes returned close to the same configuration; larger changes usually
led to worse results. No row beats the separately improved near-miss.

The bridge experiment offers individual rays rather than complete sign orbits.
Bridges, axes, and six-dimensional block sizes are all selectable. The initial
3,680-ray pool includes both exact 840 bridge constructions, integer six-block
rays, D6-type rays, and E6 spinor candidates. Its 134,000 incompatibilities are
computed exactly in `Q(sqrt(2),sqrt(3))`. Both 840 positive controls pass.
The MILP solver reported infeasibility for 841 and 842 in this finite pool;
these solver statuses are not exported formal proof certificates.

The enlarged pool has 9,424 rays and 2,187,808 exact incompatibilities. Its MILP
was terminated without a conclusion after it exceeded its requested solve
budget. Bounded native searches across the two pools made **2,489,160,192
exchanges**, without reaching compatibility at 842; the best sets still had six
conflicting pairs. Continuous refinement also tested 64 pruned 841 candidates,
without a discovery success. Finite-pool failures say nothing about unrestricted
spherical configurations.

## Evidence and replay

`campaign.json` contains completed counts and outcomes. `evidence/` retains small
certificates, coordinates, summaries, and run configurations. Large checkpoints,
graphs, and all intermediate candidates remain in ignored `runs/`.
`source_manifest.json` pins the current sources and shared dependencies;
`environment.json` records installed package versions.
The initial batch-64 wrapper is preserved in `snapshots/control_gpu_initial.py`,
and its hash agrees with the hash captured when that run started.

Every compute command was pinned to CPUs 0–126, reserving CPU 127. Concurrent
configured workers stayed below 127, with one BLAS thread per numerical worker.
GPU work used the RTX A4000. The batch-64 GPU checkpoint was deliberately stopped
and resumed on CPU after hardware benchmarking; all retained trajectories
completed their scheduled stages. No search jobs are intended to remain running.

From this directory, install the neighboring quaternionic campaign's
`requirements.txt` in your Python environment to replay the compact local
certificates using only the committed `evidence/` files:

```sh
taskset -c 0-126 python local_certificate.py
```

The following export and aggregation commands additionally require the original
ignored `runs/` data or a complete campaign rerun:

```sh
taskset -c 0-126 python certify_results.py
taskset -c 0-126 python trade_pool.py
taskset -c 0-126 python trade_pool.py --extended
taskset -c 0-126 python summarize.py
```

The standalone local certificate can be replayed from `evidence/` without the
ignored numerical checkpoints; see its proof document. Other full replays use
the neighboring two campaigns' pinned public inputs. The scripts' command-line
arguments and retained configs specify the search budgets and seeds.
