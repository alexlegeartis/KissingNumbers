# Discovery of signed P48 classes in dimensions 49–60

The frozen, authoritative result is the
[**signed-refill proof package**](../../../verifications/improved/dim49-60-p48-signed-refill/).
It supplies all twelve final bounds, the compressed certificates, complete
verification code, exact lift proof, input provenance and independent reviews.
This directory retains the discovery programs and small search records.
Reproducible binary caches, intermediate NPZ files and duplicate final witnesses
are deliberately excluded from publication by `.gitignore`.

The new mother class has **15,622 signed norm-six vectors**, including **10
without their antipodes**. Its predecessor had 15,598 signed vectors. Every
distinct pair has raw inner product at most two in the published P48p Gram
metric. This gives **52,431,624** kissing points in 49D, a gain of 24. Reusing
the exchanges, with individual compatibility and ownership checks, improves
all twelve dimensions 49–60; the new 60D bound is **58,631,571**, a gain of
3,860. These numerical gains are modest. The useful extension is that the
mother classes need not be antipodal.

## Verify the published results

From the repository root, with NumPy and a C++17 compiler:

```sh
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python verifications/improved/dim49-60-p48-signed-refill/final.py
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python verifications/improved/dim49-60-p48-signed-refill/verify.py
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python verifications/improved/dim49-60-p48-signed-refill/negative_controls.py
```

Apply the CPU-affinity requirements in the repository's `AGENTS.md` to compute
commands. `final.py` derives counts; `verify.py` proves the finite geometry and
ownership. Neither command requires running the searches below or generating
any file in this research directory.

The complete mother class was independently checked in bounded binary64 and
scalar C++ int64 arithmetic over **122,015,631 unordered pairs**. The twelve
complete final-family replays additionally checked **30,283,653** source-base
pairs and **3,540,455,163** exceptional pairs, plus isometries, exact signed
ownership, direction codes and pole codes. The published P48p minimum six and
52,416,000-point minimal shell remain explicit mathematical inputs. The full
mother shell is not enumerated.

## Reproduce discovery

The search programs need NumPy and the pinned source files already in the
repository. Run these sequentially from the repository root:

```sh
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python research/dim32plus_20261113/middim/signed_p48.py
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python research/dim32plus_20261113/middim/p48_lns.py 20261113
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/dim32plus_20261113/middim/p48_propagate.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/dim32plus_20261113/middim/compare_reproduction.py
```

The first program writes `p48_signed_class.npz`; the recorded run attained
15,614 points. The second reads that intermediate and writes
`p48_signed_class_lns.npz`, attaining 15,622. The third reads that final class
and writes the twelve `p48_signed_delta*.json.gz` files. These are generated
outputs and are not needed for publication. The final comparison checks the
complete signed point set and all twelve decoded deltas against the canonical
package, ignoring only elapsed-time metadata. It repeats no pair calculations.
Thus a successful comparison establishes that discovery reproduced the already
verified construction, rather than merely reproducing a count.

`signed_p48.py` starts with both signs of the previous 7,799-line class and
both signs of the retained low-conflict candidate pool. It makes one-for-two
signed exchanges. `p48_lns.py` searches small complete blocker neighborhoods
and finds two-for-three exchanges. The final template removes 32 signed points
and inserts 56. The finite pool and termination do not prove exhaustion.

`p48_propagate.py` applies this template separately to every old class. It
checks proposed vectors against all retained points, resolves complete signed
ownership, and accepts only positive exchanges. Different old extra points and
collisions change the attainable gains; the script does not simply multiply
24 by the number of classes. Discovery logs and small JSON histories are
retained here. The proof package's coordinate certificates are authoritative.

## Why dropping antipodality works

Let the (C_i) be disjoint signed subsets of the unit P48p minimal shell,
with within-class inner products at most (1/3). Assign a distinct antipodal
tail pair (± z_i) to each class. Each chosen signed mother point (u) is
removed from the equator and replaced by

\[
 (\sqrt3\,u/2,z_i/2),\qquad(\sqrt3\,u/2,-z_i/2).
\]

For (S=\sum_i|C_i|) and (P) poles, the count is
(52,416,000-S+2S+P=52,416,000+S+P). The pair inequalities depend on distinct
signed heads and the complete antipodal direction code. They do not require
that (-u) be in the same class, or be selected at all. The
[canonical derivation](../../../verifications/improved/dim49-60-p48-signed-refill/README.md#why-a-signed-class-is-sufficient)
gives every pair case and the ownership proof.

## Other investigated directions

`affine_profile_search.py` reconstructs the P48 minimal shell's E8 fiber
population table and searches unions of normalized affine profiles. To replay
only the four completed small-codimension experiments, with NumPy and SciPy:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/dim32plus_20261113/middim/affine_profile_search.py A1 A2 A3 D4
```

The retained JSON files explicitly mark these as candidates whose numerical
geometry needs exact replay. Their scores do not exceed the relevant
incumbents and are **not new certified bounds or maximum-clique proofs**.
The D5 run was interrupted; the larger E6–E8 cases were not completed and
have no result claimed here. Generated profile graphs are reproducible caches.

A proposed 60D route through strongly overlapping dual octads of the binary
[64,28,16] code failed its initial structural check: all 72 retained dual
octads meet pairwise in zero or one coordinate. The strongly correlated
four-normal family sought in that particular catalogue is absent. This does
not exclude general sections or other normal vectors.

## Verification runtime

The completed twelve final-family replays recorded **194.19 seconds in total**
on this session's CPU allocation. This subtotal excludes the single-class C++
check and compilation, the source-base replay and the tail/pole checks, which
were also completed separately. The complete wrapper was not separately timed.
The six-minute `run_all.py` budget is adequate for the observed workload with
one or two BLAS threads; it should not be presented as a measured six-minute
runtime. The current runner uses this value as a planning budget and supplies
no hard subprocess timeout.
