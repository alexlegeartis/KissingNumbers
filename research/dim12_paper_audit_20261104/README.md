# Audit of Takhanov–Assylbekov–Yun, arXiv:2606.18984v2

**The 841-point result survives exact verification. No 842-point arrangement
was obtained in this audit.** The useful findings are restrictions of the
structural model and the released initialization, rather than an invalid
proof of the published lower bound.

The paper and both theorem proofs were read from
[version 2](https://arxiv.org/abs/2606.18984v2), alongside the authors'
[code and data](https://github.com/k-nic/841_in_12D). Public input hashes are in
`audit.json`; the exact GitHub blob identifiers for the code are in
`source_refs.json`.

## The lower bound is genuine

`decimal841.py` reads the actual printed coordinate strings as exact finite
decimals. It proves all 353220 pair inequalities for the mathematically
normalized rays using integers. Every pairwise cosine is at most
8388607/16777216 = 1/2 - 2^-24. No spectral rank tolerance or floating point
arithmetic participates in the decision. See `decimal841_certificate.json`.

Theorem 1's norm and pairwise-inner-product calculations follow from its
stated hypotheses. Theorem 2's common rotations preserve the Gram conditions,
and its row-sum bound is sufficient for the coordinate conditions. I did not
find an error in these conditional proofs. Neither theorem claims that all
48-systems must have the displayed form.

## An exact 48-system outside the paired-floor model

The released archive has 500 numerical arrays from 160 trajectory labels.
After mathematical row normalization, **seven** satisfy both the coordinate
and pairwise conditions exactly; 493 fail at least one of them. The archive is
useful numerical seed data, but most entries should not be used as already
certified arrangements. See `archive_exact.json`.

One of the valid entries satisfies the stronger inequalities

    |v_ik| <= 1/2 - 2^-20
    v_i . v_j <= 1/2 - 2^-20    for i != j.

Its raw dyadic coordinates are kept in `strict48_counterexample.npy`, with
the integer certificate in `strict48_certificate.json`. Theorem 1's paired
floors always have sixteen pairs at cosine exactly 1/2, so this strict system
cannot be isometric to any configuration of that form. This is a counterexample
to treating that form as necessary, not a contradiction of the sufficient theorem.

Independent small motions of all 48 points remain feasible. This gives an open
240-dimensional subset of (S5)^48; it is not confined to Theorem 2's common
rotations. The proof and precise perturbation radius are in
[COUNTEREXAMPLE.md](COUNTEREXAMPLE.md).

## The released random-sign starts have only three initial geometries

For a six-coordinate Clebsch block, sign changes of the floor coordinate and
even-parity sign changes of the other five coordinates preserve the block.
Every bridge is invariant as a set under arbitrary coordinate sign changes.
Thus the only remaining data for a uniform-sign inserted point are the two
five-coordinate parities, one in each block. Exchanging the identical blocks
identifies the two mixed-parity cases.

`audit.py` enumerates all 4096 possible insertion sign vectors and finds exactly
three different full Gram multisets, with class sizes 1024, 2048, 1024. The group
action maps all starts within each class to each other. This establishes exactly
three isometry classes for the **released 841-point initialization**; it is not
a statement about all possible two-point initializations for 842.

An important caveat prevents overstating this finding. Adam is equivariant under
these signed permutations in exact arithmetic, but numerical differences can
grow rapidly. In one double-precision control, coordinate discrepancies grew
from roughly 1.6e-14 after one step to 3.9e-4 after ten steps. Consequently,
symmetry-equivalent starts can explore different numerical trajectories. It
would be incorrect to infer that all but three batches are necessarily wasted.
The issue is limited deliberate geometric diversity, not a false claim of
statistical independence.

## A different exact bridge system

The fixed five-product bridge pattern is not the only saturated 2+2 coupling.
Using all six labeled K6 factorizations gives another 720-point bridge system.
Together with the paper's exact Clebsch blocks it forms a valid 840-point code.
Its contact graph has 351136 triangles, whereas the exported start has 353696.
These exact counts prove the two arrangements are not isometric.

[BRIDGES.md](BRIDGES.md) gives the construction and proof. It also classifies
the explicitly limited family with 45 distinct 2+2 supports and all sixteen
signs per support: 4320 labeled product couplings and 720 labeled alternative
couplings, in two types under independent coordinate relabeling. This does not
cover partial sign patterns or general spherical codes. It does not refute
uniqueness of the K6 factorization up to isomorphism.

The independent placement of the two original blocks also matters. Moving the
second block's distinguished coordinate pair from (0,5) to (0,1), while fixing
the bridges, gives a valid 840-point code with a different contact-degree
distribution. `placements.json` verifies all thirty such right-block placements.

## Other code findings

The public checker accepts `maxdot <= 0.5 + 1e-7`. That predicate can accept an
invalid kissing arrangement; `check.py` supplies a spectral-reconstruction
counterexample with maximum cosine about 0.500000025, independently rejected
by an exact test. This does not undermine the published record, which has a
strict margin and passes the separate decimal certificate.

`polish_841.py` sets the actual learning rate to the printed schedule rate
divided by ten. This is a factual reporting mismatch, not evidence that its
chosen rates are mathematically wrong. The released search script also asserts
N=841, so the later 842/844 searches cannot be replayed by merely passing a
different N to that script. See `code_findings.json`.

## Tests of the exploitable restrictions

All compute was pinned to CPUs 0–126, reserving CPU 127. The GPU experiment used
the verified RTX A4000 outside the sandbox. CPU worker pools used at most 32
workers, with one BLAS thread each.

| Experiment | Starts | Outcome |
|---|---:|---|
| 49 unit vectors in the coordinate-bounded six-dimensional region, seeded from valid archive systems | 1024 | None feasible; best squared residual about 0.05147 |
| Riesz continuation with original, separately placed, rotated, and archived block seeds, then CPU polish | 128 | Best maximum cosine about 0.50782245 |
| Standard bridges, 98 free points with axes/bridges fixed first, then full release | 147 | Best maximum cosine about 0.50605622 |
| Alternative bridges, same free-point procedure | 147 | Best maximum cosine about 0.50549268 |

49 cube-bounded points plus the twelve signed axes would give 61 in R6; two
such blocks plus the 720 bridges would give 842. This conditional route was
tested, not assumed impossible. A separate 700-ray exact integer pool was
searched by MILP for 61 points. It contains an exactly checked 60-point positive
control; the solver reported no 61-point solution in that finite pool.

No run reached the 0.5 kissing threshold or improved the published 842-point
near-miss at 0.500901018602. The GPU exploration was a short float32 continuation
with double-precision CPU polishing, not a replay of the paper's much longer
float64 schedules. These negative runs do not exclude the new seed families.
`experiments.json` records the results; `runs/` keeps full local checkpoints.

## Replay and evidence

The exact construction proofs are independent of CUDA. Numerical search helpers
also use the neighboring `dim12_quaternionic_20261104` package. Install that
package's pinned `requirements.txt` for CPU work or `requirements-gpu.txt` for
GPU work, then run from this directory:

```sh
taskset -c 0-126 python fetch_inputs.py
taskset -c 0-126 python audit.py
taskset -c 0-126 python archive_exact.py
taskset -c 0-126 python strict48.py
taskset -c 0-126 python decimal841.py
taskset -c 0-126 python placements.py
taskset -c 0-126 python alternative_bridges.py
taskset -c 0-126 python couplings.py
taskset -c 0-126 python check.py
```

Large public inputs and full search checkpoints stay in ignored local folders.
Small certificates, exact-coordinate data, representative unsuccessful
candidates, and source hashes are kept alongside the code for review.
