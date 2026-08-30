# CLOSED — dimension 17: α = 192 exactly, so τ(17) = 5730 in the layered family

**This was the repository's one live lead. It is now closed, and the answer is Cohn–Li's.**

```
tau(17) = 5346 + 2*alpha,     alpha = alpha(Cay(K', W_4)) = 192     ->     tau(17) = 5730
```

```
python scripts/theta17.py       # the whole proof: ~15 s, numpy only, no solver
python scripts/frame17.py       # the frame, measured from the record (needs Cohn's data set)
```

There is no new record in dimension 17.

## What closed it

For four sections of KNOWLEDGE.md the obstruction was that `K'` was being treated as an
arbitrary 1024-vertex graph and handed to CP-SAT. **It is not an arbitrary graph.** `K'` is a linear code — closed under
XOR — so it is an elementary abelian group `F_2^10`, and `W_4` is a symmetric subset of it.
That makes `G = Cay(K', W_4)` a *translation scheme*, and Delsarte's linear-programming bound
applies verbatim. For a Cayley graph on an abelian group that bound is exactly Lovász's `θ'`,
and here it comes out at **192 on the nose**.

The certificate is six rational numbers. Index the 1024 characters `y` by two invariants — the
eigenvalue `lambda(y) = sum_{w in W_4} chi_y(w)`, and whether `y` is trivial on the fibre
subgroup `H` — and set

| (λ, trivial on H) | (30, no) | (12, yes) | (12, no) | (6, no) | (−4, no) | (−2, no) | rest |
|---|---|---|---|---|---|---|---|
| characters | 16 | 10 | 120 | 160 | 360 | 240 | 118 |
| `b_y` | 65/80 | 40/80 | 30/80 | 28/80 | 10/80 | 9/80 | 0 |
| mass | 13 | 5 | 45 | 56 | 45 | 27 | 0 |

`g(x) = sum_y b_y chi_y(x)` then has `ghat >= 0` by construction, `g(0) = 191`, and
`g(x) <= -1` at every one of the 963 points outside `W_4 u {0}` (887 of them at exactly `-1`,
60 at `-21/5`, 16 at `-9`). For any independent set `C`,

    0 <= sum_{c,c' in C} g(c-c') <= |C| g(0) - |C|(|C|-1)   =>   |C| <= 1 + g(0) = 192 .

Scaled by 80 every number in that check is an integer, so `theta17.py` verifies it with one
Walsh–Hadamard transform and no floating point at all. The matching lower bound is Cohn–Li's
own configuration — six pairwise non-adjacent fibres, 192 points, all 18 336 pairs checked.

`frame17.py` establishes that this is the right graph without using §74's derivation at all.
Sweeping all 65 536 sign patterns against the record's other 5 346 points, exactly **1024** are
admissible at height 1/3, they form one coset of a group, and that group's difference weights
are `{0:1, 4:60, 6:256, 8:390, 10:256, 12:60, 16:1}` — sixty of weight 4, which is `W_4`. It
also records something the algebra does not predict: **the record is not antipodally
symmetric.** Its 5 346-point core is, but the two height-1/3 layers are different sets, and
they have to be: one direction used at both `+1/3` and `-1/3` would sit at cos 7/9.

One wrinkle worth stating, because it looks like a discrepancy and is not. The group
`frame17.py` measures is **not the same set of words** as the `K'` `theta17.py` builds:
`theta17` builds `RM(2,4)` in the Boolean-cube coordinate order and Cohn's file uses its own,
so the two are isomorphic but differently labelled. Rather than hunt for the permutation,
`frame17.py` applies *the same six numbers* to the graph it measured. They certify it, with
the same `g`-distribution (887 at `-1`, 60 at `-21/5`, 16 at `-9`) and the same spectrum. So
the bound holds for the object the geometry actually produces, and the labelling never has to
be settled.

## Why the earlier machinery could not get there

Two natural relaxations both stop at exactly 256, and now it is clear why they are the same
number:

* the **clique cover**: `omega(Q) = 4` forces a cover of `Q` by 8 four-cliques, each of exact
  local bound 32, giving `8 * 32 = 256`;
* the **Hoffman ratio bound**: the spectrum of `G` is `60^1 30^16 12^130 6^160 (−2)^240
  (−4)^375 (−10)^96 (−20)^6`, so `n(-λ_min)/(d-λ_min) = 1024*20/80 = 256`.

Any relaxation that admits the uniform fractional point `x ≡ 1/4` on the quotient returns 256,
and clique cuts, block-pair cuts and odd-hole cuts all do. The LP over the *group* does not go
through the quotient at all, which is why it sees the extra factor.

## What this replaces

The CP-SAT route is now retired. It was a lottery, and the caveat that used to occupy this
README is worth keeping visible as a lesson rather than as a live risk — three runs of the same
capped model returned dual bounds **192, 193 and 195**, and the deterministic single-worker
variant returned **238**:

```
2400s, 8 workers                 : incumbent 138, dual bound 192
3300s, 8 workers                 : incumbent 138, dual bound 193
3300s, decision form (sum >= 193): UNKNOWN,       dual bound 195
3300s, 1 worker + fixed seed     : incumbent 135, dual bound 238
```

Each was a sound bound *for its own run*, but a headline resting on one unreproducible run that
lands exactly on the threshold is not good enough. `theta17.py` needs no such caveat: same six
numbers, same integer arithmetic, same answer on every machine. `scripts/cap22.py` and
`scripts/cap22_det.py` are kept for the record and are no longer part of the argument.

## Scope

Unchanged from KNOWLEDGE §75. This closes the **layered family** over a `Lambda_16` equator
with an antipodal pole pair and integral tiers — the family containing Cohn–Li's construction
and every usable variant of it (§74i, which showed the two extra admissible tiers are available
and worthless at a 1-to-16 exchange rate). It is not a proof that `tau(17) <= 5730` for
arbitrary configurations; the best unconditional upper bound remains the LP/SDP literature's.

What it does settle is the question this directory existed to ask: **no configuration in that
family beats 5730**, so there is no dimension-17 record to be had here.

## Files

| file | what |
|---|---|
| `scripts/theta17.py` | **the proof** — builds `K'`, checks the certificate in integers, exhibits the 192-set |
| `scripts/frame17.py` | **the frame, measured** — derives `5346 + 2α` and the graph itself from Cohn's coordinates |
| `scripts/decode17.py` | decodes Cohn–Li's record into its layers `1 + 512 + 192 + 4320 + 192 + 512 + 1` |
| `scripts/codes17.py`, `extract.py` | the tier/support enumeration behind §74 |
| `scripts/cap22.py`, `cap22_decide.py`, `cap22_det.py`, `cap28.py` | the retired CP-SAT models, kept for the record |
| `scripts/cap22*.log` | the runs those models produced — the record itself, and where the dual bounds tabulated above came from |
| `data/d17_supports.npy` | the 448 admissible 6-supports; `theta17.py` rebuilds `K'` from these |
| `data/d17_equator.npy` | the equator layer, read by `codes17.py` |
| `data/d17_tierB_best.npy`, `d17_tier_*.npy` | recorded output of the tier enumeration of §74 — 181 distinct 16-bit words in the first, the two tier levels in the others. No script here reads them; they are the record of that enumeration, kept alongside the retired models |
