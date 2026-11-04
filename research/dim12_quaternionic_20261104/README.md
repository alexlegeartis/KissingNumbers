# 12D construction campaign — 4 November 2026

**No new lower bound was obtained. The verified lower bound remains 841.**
The best 842-point maximum cosine remains the authors' existing
0.5009010186022465; feasibility requires at most 0.5. This campaign did not
improve that near-miss either. The repository's bound registry is unchanged.

The useful structural result is a **three-parameter continuous 840-point
family**, with an exact proof in [BALANCED_FAMILY.md](BALANCED_FAMILY.md).
It balances the two last quaternionic factors and motivated a search for two
simultaneous 25-point regions. This is a research construction of size 840,
not a claimed 842-point arrangement.

## Compute and experimental scope

The host exposed 128 logical AMD EPYC 7543 CPU threads, about 472 GiB RAM, and
an NVIDIA RTX A4000 with 16 GiB VRAM. CUDA was hidden inside the execution
sandbox; the GPU was subsequently verified and used outside it. The original
CPU batches used up to 128 workers, with one BLAS thread per worker. Newly
added repository guidance reserves CPU 127; all subsequent commands and replay
commands below use affinity 0–126 and at most 127 CPU workers.
A conservative displacement bound
lets the optimizer omit provably inactive pairs between neighbor-list rebuilds.
Under concurrent load, this reduced a measured objective-evaluation time from
about 29 ms to about 2 ms. It changed execution speed, not the acceptance test.

The CPU batches contain **1,568 construction trials**. Two additional GPU
batches contribute another **256 trajectories**, each polished in double
precision on the CPU. The total is **1,824**:

| Batch | Trials | Best final maximum cosine |
|---|---:|---:|
| Quaternionic 816+26, record-derived, 844-to-842, and global starts | 384 | 0.5068035031 |
| Ambient-rank lifting and decompression | 128 | 0.5009817120 |
| Alternative whole-cell 840 seeds, then deformation | 128 | 0.5013400197 |
| Cyclic 843-point orbits, prune one, release symmetry | 128 | 0.5017198170 |
| Reflected Gram feasibility | 32 | 0.5009010186, retained published seed |
| Finer eight-point-orbit candidates, prune six from 848, release symmetry | 256 | 0.5011101856 |
| Balanced family, two 25-point regions, release core | 512 | 0.5012297518 |
| GPU diameter annealing, then CPU polish | 256 | 0.5009402907 |

Every row fails the 0.5 feasibility threshold. Squared-overlap objectives and
maximum-cosine objectives differ: a lower sum of squared overlaps alone does
not mean a better kissing configuration. [campaign.json](campaign.json) and
[trial metrics](evidence/trial_metrics.json) retain both metrics.

The GPU runs use FIRE dynamics and temporary mean-zero diameter variations to
change the contact pattern. Diameters return to one for the final 40 percent
of each GPU time budget. The first run retains the smallest physical overlap
energy seen; the second carries its final contact patterns into CPU polishing.
Float32 is used only for exploration. The CUDA force was compared with an
independent double-precision CPU calculation, and the resulting candidates use
the same double-precision polish and exact verifier as the CPU batches.

The finer orbit graph has **20,736 vertices**. It uses exact arithmetic for
compatibility of eight-point independent-sign orbits, with all 72 axial points
present. Three native 128-thread tabu batches made **945,486,336 local moves**.
None found 97 mutually compatible orbits, which would have given 848 points.
The last batch saved 286 approximate sets; 256 entered the refinement above.
The full graph occupies about 411 MiB and is regenerable, rather than tracked.

An additional **8,192** searches optimized one antipodal insertion direction
and the balanced-family parameters. None found an insertable pair. The best
observed maximum cosine to the family was about 0.6035533906. This is the best
result of those searches, **not** a certified covering radius or maximality proof.

Early dense-kernel pilot batches were superseded after the faster kernel was
validated. Their partial checkpoints and stop records remain locally in
`runs/ensemble1/` and `runs/continuation1/`; they are excluded from the 1,568.

## Exact structural findings and verification

The exact whole-cell audit enumerates 1,296 possible 128-point cells and 49,248
compatible cell pairs. An exhaustive 81,643-node clique search rules out seven
cells; a six-cell witness gives 840 points with the 72 axes. This proves a
ceiling **only for the fixed-frame, fixed-weight, affine-xor, whole-cell ansatz**
specified in [exact_cells_report.json](exact_cells_report.json). It does not
exclude continuous weight changes, partial cells, or arbitrary configurations.

For the continuous family, the M2,r weights are
`(1/2, sqrt(a_r), sqrt(3/4-a_r))`, independently for `a_r` in `[1/4,1/2]`.
Its balanced frame spectrum is `78` four times and `66` eight times, compared
with `78` eight times and `54` four times at the published corner. The balance
transfers room between the two last factors. The proof checks all 2,304 mixed
unsigned label pairs using integers, and proves the remaining inequalities for
the full real parameter intervals.

`verify.py` independently checks every pair using exact Python integers.
It interprets the supplied binary64 coordinates as rays, normalized
mathematically. For a positive dot product it tests
`4*dot**2 <= squared_norm_i*squared_norm_j`, with no numerical tolerance.
It certifies the known 841-point record and rejects the published 842-point
near-miss, finding 10,937 violating pairs. Representative campaign candidates
are independently rejected in [evidence/verification.json](evidence/verification.json).

The checks in [checks.json](checks.json) include equality at cosine 1/2,
violations around 1e-13, duplicates, invalid inputs, conservative neighbor-list
coverage, independent finite-difference gradients, and clique-search comparisons
against brute force. The native tabu solver also periodically recomputes its
conflict counts from scratch.

Earlier repository notes inferred rigidity, maximality, or exact family ceilings
from contact counts and failed numerical searches. Those inferences were too
strong. This campaign supplies an exact proof for its explicitly enumerated
whole-cell family; all other negative numerical experiments remain search
evidence, rather than impossibility results.

## Provenance and replay

The initial 840 model and record alignment come from
[Takhanov's quaternionic construction](https://arxiv.org/abs/2609.09179) and its
[public coordinates](https://github.com/k-nic/quaternionic). The 842 and 844
near-misses come from the
[authors' earlier repository](https://github.com/k-nic/841_in_12D/tree/main/beyond841).
[inputs.json](inputs.json) pins SHA256 hashes and baseline measurements.
[source_manifest.json](source_manifest.json) pins the final search sources.

From this directory, using Python 3.13 and the pinned dependencies:

```sh
python -m pip install -r requirements.txt
taskset -c 0-126 python fetch_inputs.py --include-near-misses
OPENBLAS_NUM_THREADS=1 taskset -c 0-126 python prepare.py
taskset -c 0-126 python check.py
taskset -c 0-126 python exact_cells.py
taskset -c 0-126 python verify.py data/record841.npy
```

The archived batch configurations in `evidence/*_config.json` specify the
trial counts, seeds, step budgets, and source hashes. For example:

```sh
taskset -c 0-126 python balanced_family.py --trials 512 --workers 112 --steps 4000 --rounds 3 --out runs/balanced_family
taskset -c 0-126 python family_holes.py --trials 8192 --workers 16 --out runs/family_holes
taskset -c 0-126 python quaternion_cells.py
taskset -c 0-126 python orbit_pool.py --threads 96
taskset -c 0-126 python orbit_weights.py
taskset -c 0-126 g++ -O3 -march=native -std=c++17 -fopenmp orbit_tabu.cpp -o runs/orbit_tabu_partial_bin
taskset -c 0-126 runs/orbit_tabu_partial_bin runs/orbit_pool runs/replayed_orbit_tabu 180 127
```

Native runs use wall-clock budgets, so iteration counts and saved-file ordering
depend on scheduling. The saved candidate point sets and their hashes are the
reviewable artifacts. `runs/` keeps the complete local logs, larger inputs,
graphs, and checkpoints; representative arrays and compact reports are kept
in `evidence/` as small, versionable artifacts.

GPU replay additionally needs `requirements-gpu.txt`, CUDA device access, and
`CUPY_CACHE_DIR` pointing inside the workspace. Select the RTX A4000 at device
0 as verified on this host. Run `gpu_fire.py` with CPU affinity 0–126; the script
refuses this host's unreserved CPU affinity and an unexpected GPU model.

The resulting evidence favors searching for a different mixed-core geometry or
coupling between the two exceptional regions. Increasing the polishing budget
on the same published 842 contact pattern has not demonstrated progress here.
