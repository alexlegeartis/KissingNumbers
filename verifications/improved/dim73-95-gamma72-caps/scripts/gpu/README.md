# Gamma_72 classes: a standalone GPU package

Plain `numpy` + `torch`. No Colab, no notebook, nothing device-specific -- the exact-arithmetic
back end is chosen by probing the device and checking it against an int64 reference, and there
is a float32 fallback that runs on any GPU and on the CPU.

```bash
pip install -r requirements.txt
python run.py --auto --classes 31200 --seed 7 --tag prod          # ~2 h on a T4
python run.py --auto --classes 31200 --seed 7 --tag prod --resume # after any interruption
python checkpool.py out/cert_prod.json                            # the pool has no repeated line
python regenerate.py out/cert_prod.json                           # rebuild from the seed, match the hash
```

`--auto` reads the free memory of whichever device it finds and sizes the pool and the scan
window to it, keeping the total under 60% (see the memory note below -- this is the one tuning
that matters). Everything else has a working default. Output is three files in `--out`
(default `./out`): a ~90 kB certificate, a ~120 kB array of class sizes, and an ~11 MB
resumable state.

A run can be interrupted at any point and continued with `--resume`: the state is an alive-line
bitmask plus one SHA-256 digest a class, every class is verified exactly as it is built, and
disjointness is a running invariant (`|alive| + sum(sizes) = |pool|`). A resumed run reproduces
an uninterrupted one bit for bit -- same class sizes, same final hash -- **on the same software
stack**; across stacks it does not, and the state now carries a pool digest that a resume checks
before doing anything (see *Reproducibility, and its limit* below).

## The quantity

Dimension `72+k`'s cap construction over Gamma_72 gains

    gain(k) = 4 * S(T) + 2 * (S(T+P) - S(T)),      S(j) = sum of the j largest class sizes

over a family of **pairwise-disjoint classes**, a class being a set of Gamma_72 minimal LINES
with pairwise `|<u,v>| <= 2`.  `T` and `P` come from the zero-sum triple partition of a kissing
configuration of `R^k`; for `k = 23` they are `30 990` and `82`.  So the whole question is how
many minimal lines fit into ~31 072 disjoint classes, and nothing else about them matters.

Greedy on a candidate set of `n` lines returns about `ln(1 + n p)/p`, where `p` is the density
of conflicting pairs.  Both inputs are therefore worth exactly what their logarithm says:

| pool                                | n      | p        | class |
|-------------------------------------|--------|----------|-------|
| one automorphism orbit               | 2.06M  | 0.010241 |   972 |
| orbit + difference rule (CPU, 2026)  | 3.39M  | 0.005469 |  1797 |
| this module                          | 90M    | 0.004802 |  2175 |

and `p` has an exact floor: the minimal vectors of Gamma_72 are a spherical 11-design, which
pins the inner-product distribution down completely,

    <u,v> = 0 : 2 603 658 750    +-1 : 1 512 243 200 each    +-2 : 280 928 256 each
           +-3 :    13 959 168 each    +-4 :     127 800 each    +-8 : 1

so a uniform sample of the shell has `p = 28 173 936 / 6 218 175 600 = 0.0045309`.  `run.py`
prints the measured histogram against this exact one; it is the pool's quality certificate.

## Two moves, and each does what the other cannot

* **walk** -- push a random sample of the pool by an independent random word in the
  generators.  Uniform measure is stationary for a random walk on a finite group, so a long
  word lands close to uniformly in each line's own orbit.  This FILLS orbits and can never
  leave one; `|Aut|` caps an orbit at about 2.06M of the 3 109 087 800 lines, so the walk
  alone saturates -- a first version stalled at ~25M.
* **diff** -- the production rule `<u,v> = +-4  =>  u -+ v is minimal`, swept over the whole
  pool.  This is the only move that OPENS a new orbit.  Sweeping one `u` against the whole
  pool finds all `n * P(|ip| = 4)` of its neighbours at once, where sampling `k` lines and
  taking the pairs inside the sample finds only `C(k,2) * P`; the difference is between a
  sweep that returns 3.9M lines and one that dries up.

Alternating them reaches 90M lines in 87 s and is still growing (+3.9M a sweep at the end).
The walk sample must be SMALL: the pool is a union of isometric copies of the sample, so the
within-copy pairs carry the sample's own density with weight `1/#copies`.

## Exactness

`<x,y> = (x G) . y`.  Gamma_72 is unimodular, so `G^-1` is integral and, for minimal `x`,

    |x_j|     = |<x, b*_j>| <= sqrt(8 * max_j (G^-1)_jj) = sqrt(8 * 72) = 24 ,
    |(x G)_j| = |<x, b_j>|  <= sqrt(8 * max_j G_jj)      = sqrt(8 *  8) =  8 ,

hence `|<x,y>| <= 72 * 24 * 8 = 13 824`.  `check_bounds()` recomputes both from the Gram at
start-up so the constants cannot drift -- an earlier version assumed 16 and 4 from a sample of
the working pool and the walk produced a coordinate of 20 within seconds.

Three back ends are tried and each is checked against an int64 reference before any class is
built: `torch._int_mm` (int8 tensor cores, exact by definition), float16 in with float32
accumulation (inputs and the only test made, `|ip| <= 2`, are exact; anything of modulus
>= 2048 comes back within 8, so it is never mistaken for one of modulus <= 2), and float32
(exact to 2^24).  **TF32 must be off** -- it keeps 10 mantissa bits, exact only to 1024.

## The certificate

Every random decision is made by numpy's PCG64 from `--seed` and only then moved to the
device, which does nothing but exact integer GEMM.  The run is therefore reproducible from the
seed **within one software stack** -- see *Reproducibility, and its limit* below for what that
does and does not cover -- and the 6.5 GB of lines never has to be shipped:

    cert_<tag>.json    parameters, seed, |ip| histogram, class sizes, SHA-256, timings  (~90 kB)
    sizes_<tag>.npy    class sizes, int32                                              (~120 kB)

`regenerate.py cert_<tag>.json [n]` rebuilds from the seed, re-verifies in exact integer
arithmetic and compares the hash; the optional `n` checks a prefix of `n` classes only.



## What a bigger card buys

`--auto` sizes both the pool and the window, and the class grows as `ln(1 + scan*p)/p`, so the
return is logarithmic in memory but it is real -- and dimension 95 reads 4x the sum over 30 990
classes, so a hundred lines a class is 12M on the bound:

| free memory | pool | scan | class (predicted) | dimension 95 |
|---|---|---|---|---|
|  8 GB |  44M | 4.4M | 1895 | 6 453 million |
| 15 GB (T4) |  87M | 8.7M | 2025 | 6 469 million |
| 24 GB | 142M | 14.2M | 2119 | 6 481 million |
| 40 GB | 241M | 24.1M | 2219 | 6 493 million |
| 80 GB | 400M | 40.0M | 2315 | 6 505 million |

The T4 run these were calibrated against used pool 80M, scan 8M and realised a mean of 1997
against the 2025 predicted, so the column is about 1.5% optimistic.  Time scales with `--scan`
too (0.23 s a class at 8M on a T4), so a 40 GB card is roughly 3x the wall clock per class as
well as +10% on the class -- worth it if the card is otherwise idle.

## The memory rule (the one tuning that matters)

A class holds the pool (int8), the float16 window, the copy the compaction makes of it, the
alive index and the product tile.  Measured on a T4 with 14.7 GB:

| residency | seconds a class |
|---|---|
| 9.2 GB (pool 80M, scan 8M)  | **0.23** |
| 11.5 GB (pool 88M, scan 16M) | **6** |

Nothing about the arithmetic changes between those two rows; past roughly 70% of what is free
the caching allocator starts evicting and re-serving the same buffers.  `--auto` targets 60%.
It took most of a session to find, because it looks exactly like a slow GPU -- a second card
reproduced it, which is what ruled the hardware out.  `run.py` prints the projected residency
and warns if it is above 72% of free memory.  If a run is unexpectedly slow, lower `--scan`
first; the class shrinks only as `ln(1 + scan*p)/p`, so halving it costs about 6%.

## Measured on the T4

| what                                  | rate                                       |
|---------------------------------------|--------------------------------------------|
| pool, 88M lines (seeds + diff + walk) | 85 s                                       |
| classes, scan 10M, pool 90M           | 0.39 s each, sizes 2015..2079 (mean 2047)  |
| classes, scan 16M, pool 88M           | **6 s** each -- see below                  |
| exact re-verification                 | 5 ms a class (int8 tensor cores)           |

**Stay under about 10 GB.**  This is the one tuning that matters and the only one the model
missed.  A class at scan `s` holds the pool (80M x 72 int8 = 5.8 GB), the gathered window
(`s` x 144 in float16), the copy the compaction makes of it, and the product tile.  At scan
10M that is 9.2 GB and a class takes 0.39 s; at scan 16M it is 11.5 GB and a class takes 6 s.
Nothing about the arithmetic changed -- the caching allocator starts evicting and re-serving
the same buffers.  It cost most of a session to find, because the symptom looks like a slow
GPU: the same run on a second T4 was equally slow, which is what ruled the hardware out.

Three other tunings that were tested and did **not** matter:

* tiling the `(alive x block)` product into the T4's 4 MB L2 (`pchunk`) is 2.8x SLOWER at 4096
  than at 2M -- launch overhead beats the traffic saved, so the product is not L2-resident the
  way the model assumed;
* whether the alive list is left in pool order or randomly permuted: 1920.4 against 1920.2
  lines a class over the first 1000, on the same pool and seed.  Identical.  (Two earlier
  readings that seemed to show otherwise were both wrong: one A/B toggled a flag that an
  edit had dropped, so both arms permuted, and a later comparison was against a run using a
  LARGER scan window, not a different ordering.  What actually differed was `--scan`.);
* an int8 operand halves the gather and compaction traffic, but `torch._int_mm` returns int32,
  which doubles the product traffic -- and the product dominates.  float16 in, float32
  accumulate is both exact and 2x faster (0.76 s a class against 1.57 s at scan 20M).

`torch._int_mm` also insists on inner and output dimensions that are multiples of 8 and at
least 32 rows, which the accepted block of a greedy round does not respect (73 on the run that
caught it); `mm()` pads with zeros and slices.  The CPU implementation does not enforce this,
so it only ever shows up on the device.

## Reproducibility, and its limit

Every random DECISION is numpy's PCG64 from `--seed`, and the device does nothing but exact
integer GEMM -- which is what lets the certificate be 90 kB instead of the 6.5 GB the lines
would take.  It is **not** the same as reproducibility on any machine.  The pool is
*assembled* by torch primitives whose tie-breaking is unspecified (`dedup_into` sorts by a
64-bit key), and the pool's ROW ORDER is what the greedy window is drawn from, so two software
stacks can build two different pools from one seed.  Measured 2026-08-25: a certificate made
on a T4 under one Colab image did not regenerate under a later one -- the first 300 classes
came back summing to 595 325 against the certificate's 601 370 -- while two builds inside one
process agreed bit for bit.  `diag_repro.py` separates the four hypotheses and settles three
of them (the pool IS reproducible in-process, `ip_histogram` is inert, writing the state is
inert), which leaves the environment.

**This was dangerous, not untidy.**  Everything downstream is indexed by pool POSITION: the
alive bitmask a `--resume` restores, the per-class digests, the family SHA-256.  A resume
against a different pool silently reuses lines already taken, and the running invariant
`|alive| + sum(sizes) = |pool|` counts rows, so it would not notice.

So `gamma72.pool_digest()` fingerprints the pool -- SHA-256 of every 79th row, in order --
`run.py` records it in the certificate, `build_classes` stores it in the resumable state and
**refuses to resume** without a match, `checkpool.py` asserts it and `regenerate.py` reports
it.  A state written before the guard existed carries no digest and is refused rather than
trusted.  Verified end to end by resuming a run from a fetched state on a fresh VM: the digest
matched and it picked up at exactly the right class.
