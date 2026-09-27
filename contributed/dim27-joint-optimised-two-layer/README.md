# Dimension 27: the two-layer configuration, jointly re-optimised

**τ(27) ≥ 201 566**, against the **201 509** of [`verifications/improved/dim26-27-iota-triangles`](../../verifications/improved/dim26-27-iota-triangles/)
at commit `c349d56`. An improvement of **+57**, by optimisation alone: the construction, and every
mechanism in it, is that package's.

```
pip install -r requirements.txt
python verify27_param.py 201256                                                          # ~2 min
PYTHONPATH=independent python independent/verify2627_layers_independent.py . 27 201566     # ~2 min
```

Expected last lines:

```
ALL CHECKS PASS   K(27) >= 201566
ALL CHECKS PASS      K(27) >= 201566   (independent, exact, two layers)
```

Both read only `data/` and rebuild the 196 560 Leech minimal vectors from the Golay code, each from
its own implementation. They are two separately written checks: the first is this repository's own
`verify27.py`, the second was written from the geometric description and shares no code with it.

## The configuration

The format is exactly that of `dim26-27-iota-triangles/data/`; the four files can be copied over the
ones there.

| | triangle 0 | 1 | 2 | 3 | total |
|---|---|---|---|---|---|
| class heads `y = 3u + v` (6 leans) | 579 | 437 | 622 | 452 | **2090** |
| free heads `y = ±2v` | 39 | 45 | 41 | 43 | **168** |

| | line 0 | line 1 | line 2 | total |
|---|---|---|---|---|
| second-layer heads `(√3/2)u` at `|x|² = 3` | 105 | 104 | 101 | **310** |

| family | count |
|---|---|
| equator `(z, 0)`, `z` a minimal vector that is not an owner | 194 160 |
| caps, three per first-layer head | 3 × 2258 |
| caps, two per second-layer head | 2 × 310 |
| axis | 12 |
| **total** | **201 566** |

    count = 196560 - 2090 - 310 + 3 * 2258 + 2 * 310 + 12 = 201566          (the first layer alone gives 201256)

The shipped configuration has 2133 + 129 first-layer heads and 284 second-layer heads. This one gives up
43 class heads for 39 more free heads and 26 more second-layer heads.

## How it was found

Every item type is indexed by minimal vectors. A class head is `y = 2u + w` for an ordered pair `(u, w)`
of minimal vectors with `⟨u, w⟩ = 1` (norm-4 units) — 196 560 × 47 104 = 9 258 762 240 pairs, which is
every class head over every norm-6 lean `v = w − u`; a free head is `±2v`; a second-layer head is an
owner `u` on one of the three lines. For a head `y'` already present, `⟨y, y'⟩ = 2⟨u, y'⟩ + ⟨w, y'⟩`, so
one integer matrix (minimal vectors × present heads) decides every conflict. A GPU kernel scans all
9.26 × 10⁹ pairs against a configuration in about 3.5 minutes in exact integer arithmetic, and every
candidate it keeps is re-derived on the CPU in 64-bit integers.

A joint optimiser then lets the three item types move **together**: candidates with few blockers, the
blockers themselves, and the pairwise conflicts go into one exact 0–1 programme with weights 2 (class
head), 3 (free head) and 1 (second-layer head), solved by CP-SAT under a time limit; the scan is
repeated after every gain. The layers have to be optimised jointly: on a first layer optimised alone
(2135 class + 136 free heads, 201 250) only **105** of the 284 shipped second-layer heads still fit,
and re-packing the second layer over that fixed first layer stops at 234. Started from the shipped
201 509, from that first layer, and from the hybrid, the joint search reached 201 533, 201 545 and
201 540; recombining them and continuing gave this configuration. None of these runs proves
optimality, and the search had not converged when this was cut.

## What has been checked, and what has not

Checked, exactly, by both programmes: every head removes at most one minimal vector and the owners are
distinct; every first-layer head pair against the 48 / 96 thresholds; every head against the axis; the
axis against itself; every (first-layer head, second-layer owner) pair against `⟨Y, u⟩ ≤ 32`; every
second-layer pair (`⟨u, u'⟩ ≤ 8` on one line, `≤ 16` across lines); second-layer owners minimal,
distinct and disjoint from the first layer's; the count. The independent programme also derives the
bar 32 as `⌊32√3 − 32/√2⌋`, checks that free heads have the form `2v` with `v` of norm 6, the
distinctness of all points and the full ambient rank, and it rejects fourteen deliberately corrupted
artefacts (owner too close to a first-layer head, two heads at `⟨u,u'⟩ = 2` on one line, a shared
owner, a duplicated owner, a non-lattice owner, a bad line index, a head moved to another triangle).

**Not checked:** optimality; and **no second party has yet verified this configuration** — both
programmes were run by the contributor only.

## Files

    verify27_param.py            a copy of ../../verifications/improved/dim26-27-iota-triangles/verify27.py with
                                 one line changed (stated at its top): the pinned first-layer total is an argument
    lib/                         golay.py, leech.py, layered.py: byte-identical copies of that package's lib/
    independent/                 verify2627_layers_independent.py, verify2627_independent.py (which it imports),
                                 kiss_ref/ (a separate Golay / Leech generator)
    requirements.txt             numpy, sympy, mpmath, pinned
    data/heads27_Y.npy           the 2258 first-layer heads (int64, Cohn units, norm 192)
    data/heads27_side.npy        their triangle (0 to 3)
    data/heads27_layer2_u.npy    the 310 second-layer owners (int64, norm 32)
    data/heads27_layer2_line.npy their line (0 to 2)

SHA-256 of the data:

    e137a6e3b11e6f70a0bc631b474f4ea340078f453893c56d5452513da5264bd1  data/heads27_Y.npy
    b0e1457e2793db1ddd2fcc23a9110f8f956757ee08e017fac99200cae5199fb2  data/heads27_side.npy
    8f722bf903f5324e3e4f35f9c271d62949d39d4e56a21b0ca58b6c93bd7a54fa  data/heads27_layer2_u.npy
    1aa7054bf9f85ce9b52c887481607348dcd117016d991ab56148b869921ebf4f  data/heads27_layer2_line.npy

## Provenance

The construction is that of `dim26-27-iota-triangles` — the coset triangle on every triangle of
directions, the free-head layer and the second cap layer are all A. Kravatskiy's (joint work in
progress with H. Cohn and B. Lindow). The contribution here is the optimiser only. Contributed by
B. Lindow; computations carried out with Claude Code under his direction.
