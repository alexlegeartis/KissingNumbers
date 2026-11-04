# Larger-compute 19D and 21D campaign, 4 November 2026

No improved lower bound was found. The verified comparison configurations are
**12460 in 19D** and **31560 in 21D**, both due to Kevin Li's
[14 October paper](https://hexagonmath.org/pdf/2610.00248v1).
These supersede the older cached Ho and Cohn–Li values. The still older counts
on the accompanying GitHub repository were also excluded as benchmarks.
See `source/ATTRIBUTION.md`.

The aggregate is `data/summary.json`; `pending_runs` must be empty before the
run is considered complete. Numerical candidates are explicitly invalid until
an exact certificate proves otherwise. No main result registry was changed.

## Verified starting configurations

`prepare19.py` rebuilds the 10668-point section base, 768 marked axis vectors
and 1024 diagonal vectors from the paper's generators and Table 5. The base
and axis rays are integral; the diagonal rays have coordinates in
`Z[sqrt(6)]`. All compatibility decisions are exact. Large intermediate
squares use arbitrary-precision integers when int64 would overflow. The
reconstruction checks 135348112 matrix entries including symmetric duplicates,
with zero violations.

`prepare21.py` rebuilds the punctured Golay code, its 21 weight-five projective
lines, the published linear labeling, the 27720-point base and 3840 added
vectors with magnitudes two and three. It verifies the geometric labeling
rule on all 1806336 codeword/profile combinations. `verify21.py` independently
checks all **498001020** unordered pairs of the 31560 integer rays, with no
symmetry reduction and no violations. Small integer matrix products computed
by float64 BLAS are exact; all subsequent angle decisions use integer squares.

## 21D searches

The published line-weighted construction contains 3840 added points. The
paper's counting bound permits up to 3908, so a search should optimize labels
beyond the published linear-coset choice.

The weight-five differences span a 1024-word component. The code has four
isomorphic components, and points from different components are compatible.
For a selected word `w` with label `L`, its neighbor `w+L` must be omitted.
Along any selected word pair at difference `P`, the two labels and `P` must
be three distinct lines. These conditions are equivalent to the geometric
integer inner-product check.

On one component this is an independent-set problem with **21504** candidates
`(word,line)` and degree **81**. The current solution selects 960. Selecting
961 would give 31564 points after replication across all four components.
The graph contains every label, with no orbit-invariance restriction.

Searches completed:

| Method | Scope | Outcome |
| --- | --- | --- |
| Parallel ruin/recreate | Full component | 66159946 attempts, best 960 |
| Ruin/recreate plus exact one-for-two augmentations | Full component | 7640717 attempts, best 960 |
| CP-SAT | Full nonlinear labeling, 22528 Boolean variables, 600 seconds and 24 workers | Kept 960; returned FEASIBLE, not OPTIMAL. Its bound was 49 omitted words, leaving a gap |
| CP-SAT decision search | Two sampled sign subgroups, ranks three and two | Reported INFEASIBLE for a gain under each imposed symmetry; this does not rule out asymmetric labelings |

`solve21.py` fixes one omitted word at zero by translation, which loses no
generality. In the symmetry models, the subgroup is sampled inside the
published kernel; it is a genuine additional restriction. Solver infeasibility
was not exported as a separately replayable formal proof.

The A4000 search fixes the full base and moves the 3840 auxiliary vectors
plus one inserted point. This is more efficient than repeatedly computing
base-base interactions. It explored 23 starts, followed by a continuation
using quadratic and sixth-power penalties. The best saved candidate still
exceeds the required maximum inner product 1/2. Exact final values and
independent CPU checks are in `data/summary.json`.

## 19D searches

The published construction already takes one point for each of its 1792
admissible projected heads, so selecting different words inside that same
shape family cannot improve the count. Two broader searches were run.

First, 2048 LP-ascent starts sampled vertices of the polar of the complete
10668-point base. Rational reconstruction, exact admission checks, and the
order-16 sign action produced **3568 distinct additional rational directions**.
Together with the published auxiliary layer, this gives a 5360-point pool.
All pool-pair decisions were rebuilt exactly in `Z[sqrt(6)]`; the resulting
graph has 215072 undirected edges. The initial conservative numerical graph
and the rebuilt exact graph agree.

A parallel search with one-for-two augmentations completed **9526156** moves
without exceeding 1792 auxiliary points. CP-SAT reported INFEASIBLE for 1793
on the exact pool. This closes that sampled finite pool, not the full
continuous region. A separate attempt to extract a simple 1792-clique
partition did not succeed; no such independent clique certificate is claimed.

Second, the GPU moved all 1792 added vectors and one inserted point while
keeping the base fixed. It tried 41 starts, including reflected tails and
LP insertion directions. None produced a valid 12461-point configuration.

## Replaying

Commands assume this directory as the working directory. Dependencies are
NumPy, SciPy, OR-Tools 9.15 and PyTorch 2.6, installed in the Python environment
used below. The original campaign used OR-Tools from the enclosing workspace's
`.compute-deps`. All compute commands reserve
logical CPU 127 as required by the repository's `AGENTS.md`.

```sh
export OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=1
taskset -c 0-126 python prepare19.py
taskset -c 0-126 python prepare21.py
taskset -c 0-126 python verify21.py
taskset -c 0-126 python solve21.py --rank 0 --seconds 600 --workers 24 --seed 21
taskset -c 0-126 python solve21.py --rank 3 --seconds 300 --workers 8 --seed 2123 --gain
taskset -c 0-126 python solve21.py --rank 2 --seconds 300 --workers 8 --seed 2122 --gain
taskset -c 0-126 python polar19.py --build --exact
taskset -c 0-126 python solve_pool19.py --exact --seconds 180 --workers 16
taskset -c 0-126 python verify_results.py
```

`polar19_rays.npy` is the saved exact input for the polar-pool replay. The
random numerical sampling is not needed to verify those directions or their
exact conflict graph. `search21.cpp` is the standalone augmented local-search
engine. GPU scripts use float64 and independently tested analytic gradients;
they require supported escalated execution because CUDA is hidden inside
the sandbox. Recheck device identity before launching; only the RTX A4000
was used, and no GPU job overlapped another in this campaign.

The manifest hashes the inputs and source. Regenerable graphs, executables
and uncertified numerical candidate arrays are ignored by `.gitignore`.
