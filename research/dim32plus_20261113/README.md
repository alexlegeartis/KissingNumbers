# Kissing-number search in dimensions 32 and above — 13 November 2026

This campaign produced **thirteen verified lower-bound improvements**. The
substantial result is **tau(69) >= 743,522,304**, an increase of **54,069,382
points (7.8424%)** over the preceding repository certificate. Signed P48
classes give twelve smaller improvements in dimensions 49–60.

The comparisons below use the immediately preceding repository values, rather
than older published baselines. The dated source review is in
[theory/source_audit.json](theory/source_audit.json). The bounded literature
review does not establish exhaustive worldwide priority.

## Verified results

| Dimension | Previous repository bound | New verified lower bound | Gain |
|---|---:|---:|---:|
| 49 | 52,431,600 | **52,431,624** | +24 |
| 50 | 52,462,800 | **52,462,872** | +72 |
| 51 | 52,509,596 | **52,509,740** | +144 |
| 52 | 52,603,116 | **52,603,402** | +286 |
| 53 | 52,727,584 | **52,728,054** | +470 |
| 54 | 52,975,804 | **52,976,623** | +819 |
| 55 | 53,392,162 | **53,393,546** | +1,384 |
| 56 | 54,260,992 | **54,263,398** | +2,406 |
| 57 | 54,757,822 | **54,760,693** | +2,871 |
| 58 | 56,266,316 | **56,270,076** | +3,760 |
| 59 | 56,948,468 | **56,952,414** | +3,946 |
| 60 | 58,627,711 | **58,631,571** | +3,860 |
| 69 | 689,452,922 | **743,522,304** | +54,069,382 |

### 69D: compress four normal coordinates into one

The [Gamma72 head-compression package](../../verifications/improved/dim69-gamma72-head-compression/)
selects ten profile families of the minimal shell of Gamma72. It projects each
family onto a common 68-dimensional complement, rescales that block, and adds
one rational height. Exact angle inequalities permit five populous profiles
of squared projection norm one quarter to coexist with five further families.

Two standard-library verifiers independently reconstruct the population dual
and all 55 geometric pair types, using integral Gram coordinates and orthonormal
coordinates respectively. A third mathematical/arithmetic review and semantic
corruption controls also pass. Orientation averaging establishes existence of
a favorable tuple of actual lattice normals; the proof does not enumerate the
hundreds of millions of image coordinates or assume ambient shell symmetry.

### 49–60D: allow signed mother classes

The [signed P48 refill package](../../verifications/improved/dim49-60-p48-signed-refill/)
replaces the antipodal-class restriction by signed classes with raw inner
products at most two in squared-norm-six coordinates. The new mother class has
15,622 vectors, including ten without their antipodes. Two independent exact
arithmetic implementations check all 122,015,631 unordered pairs.

Propagation to each dimension verifies the actual final classes, distinct
signed ownership, complete antipodal direction codes, and separate pole codes.
The count is 52,416,000 plus the total number of selected signed heads plus
the pole count. The gains are modest; they are not obtained by multiplying
the single-class gain without checking the other classes.

Both constructions retain the explicitly cited published facts about their
parent lattices. They prove lower bounds for unrestricted kissing numbers;
they make no lattice-packing or global-optimality claim.

## Reproduce the final proofs

Run from the repository root, with CPU affinity appropriate to the host and
at least one available logical CPU left free, as required by `AGENTS.md`.
The 69D commands need only Python 3's standard library. The signed P48
verification additionally needs NumPy and a C++17 compiler.

```sh
python verifications/improved/dim69-gamma72-head-compression/verify.py
python verifications/improved/dim69-gamma72-head-compression/verify_orthonormal.py
python verifications/improved/dim69-gamma72-head-compression/rejection_controls.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verifications/improved/dim49-60-p48-signed-refill/verify.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verifications/improved/dim49-60-p48-signed-refill/negative_controls.py
```

Use the final packages for proof replay and registry values. Files under this
research directory preserve discovery, numerical attempts, and review history.
**Raw generated candidate pools, solver logs, floating-point duals, and numerical
search scores are research outputs, not certificates.** A saved optimizer
status is not a substitute for the exact checkers above.

## Search coverage and limits

| Lane | Work completed | Outcome |
|---|---|---|
| [lowdim](lowdim/README.md) | Large simultaneous support replacements in 34–35D; a complete specified sparse-head screen in 38D; independent 69D review | No new 32–39D bound; 69D review passed |
| [middim](middim/README.md) | Signed P48 exchanges and propagation; affine P48 profile searches; a dual-octad structural screen | Twelve verified 49–60D improvements; no new 40–48D bound |
| [gamma](gamma/README.md) | D4 and D5 affine profile models; independent orthonormal 69D certificate; finer one-head search | Recovered the existing 68D affine value and independently proved the new 69D value |
| [theory](theory/README.md) | New head-compression construction; one- and two-coordinate head searches; 51 bounded continuous 70D attempts; source audit | Certified the 69D gain; neighboring-dimension searches did not improve their records |

The experiments cover explicit finite pools, selected normal sublattices, and
bounded numerical searches. Unsuccessful runs do not prove that those entire
construction families are exhausted. The campaign is frozen with the final
certificates and retained negative search evidence above.
