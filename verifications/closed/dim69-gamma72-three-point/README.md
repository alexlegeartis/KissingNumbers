# SUPERSEDED for dimension 69 — see [`../../improved/dim68-69-gamma72-cross-sections/`](../../improved/dim68-69-gamma72-cross-sections/)

> **This page's conclusion about dimension 69 was wrong, and is kept only as a record.**
> It argued that 527 597 644 must not be claimed because the Gram (4,−2,0) could not be
> proved realizable. Realizability was being sought *through the LP's own witness lemma*,
> but it is a fact about the lattice, and Γ₇₂ is available in coordinates — the triple is
> simply **exhibited** (`realize_k3.py`), as is a D₄ four-tuple that gives dimension 68.
> **τ(69) ≥ 627 822 180 and τ(68) ≥ 361 275 480 are both claimed.**
>
> What survives here is the analysis of dimensions **53–55**, and the record of how a
> sufficient condition for realizability was mistaken for a necessary one.

Dimensions 70 and 71 are large improvements
([`../../improved/dim70-71-gamma72-cross-sections/`](../../improved/dim70-71-gamma72-cross-sections/)).

## Dimension 69: the LP says 527 597 644, and that number is not usable

The floor in dimension 69 is 331 737 984 (dimension 64, Edel–Rains–Sloane, inherited by
monotonicity). A three-point cross-section of Γ₇₂ would beat it comfortably — *if* the Gram
matrix were realizable.

| Gram (off-diagonals) | det | LP certifies | realizable? |
|---|---|---|---|
| (0, 0, 0) — three pairwise orthogonal | 512 | **149 415 783** | **YES** |
| (1, 0, 0) | 504 | 245 233 175 | not proved |
| (2, 0, 0) | 480 | 346 628 061 | not proved |
| (3, 0, 0) | 440 | 409 062 220 | not proved |
| (4, 0, 0) | 384 | 469 123 586 | not proved |
| (4, −1, 0) | 376 | 500 027 708 | not proved |
| (4, −2, 0) | 352 | **527 597 644** | **not proved** |

> **The best Gram whose realizability can be *proved* certifies only 149 415 783 — a factor
> 2.2 below the floor. So dimension 69 is not a record.**

Realizability is proved the way it is everywhere else in this project (Lemma 4 of
[`common/PROOF-kpoint.md`](../../../common/PROOF-kpoint.md)): the *k*-tuple exists if the
witnessing cell of the (k−1)-point system has strictly positive LP minimum. For every Gram
above except the orthogonal one, that cell has LP minimum **0**, so nothing is witnessed.

**This is the one place in the project where a large, tempting number had to be thrown
away.** An early results table listed 527 597 644 for dimension 69; it should be
disregarded. The prose of that same table already said dimension 69 does not reach the
record — the row was simply stale.

## Why the three-point LP is weak here

The k = 3 LP loses a factor 7 relative to k = 2 (1.07 × 10⁹ → 1.5 × 10⁸), whereas on the
Leech it loses only 2.2 (43 164 → 19 530 against a truth of 19 962). So most of the gap is
**LP slack, not geometry** — the bound is bad, not the configuration.

**Pair-marginal constraints do not help.** Feeding the 649 two-point cells of Γ₇₂ (LP minimum
and maximum for every cell of every pair inner product) into the three-point LP as bracketed
marginals — encoded as equalities with bounded slack variables, so the exact certificate
machinery still applies — changes nothing. Validated on the Leech: k = 3 orthogonal stays at
19 530 with or without them. **The two-point brackets are implied by what the three-point LP
already knows.**

## A bug worth recording

An earlier sweep keyed Gram matrices by the *multiset of |entries|*, which wrongly identifies
sign-inequivalent Grams — they have different determinants and different counts. Fixed by
enumerating realizable Grams directly from the two-point cell table.

## Dimensions 53–55 are dead for a different reason

They would come from the extremal 56-dimensional lattice, whose cross-section sizes are

```
dim 53   7 652 880
dim 54   3 888 000
dim 55   1 820 399
```

all far below the **52 416 000** that dimensions 49–63 inherit from dimension 48. So the
cross-section route is simply the wrong tool there; the cap construction over P₄₈ is what
wins, and it does — see
[`../../improved/dim49-63-p48-caps/`](../../improved/dim49-63-p48-caps/).

## Dimensions 65–68: a cheap partial win, recorded but not claimed

Dimensions 65–69 inherit 331 737 984 from dimension 64, and the Γ₇₂ cross-sections fall below
that by k = 3. The cap construction is the only lever, and it applies to **any** antipodal
kissing configuration on the equator — including the Edel–Rains–Sloane dimension-64
configuration itself. That configuration's four families are, with |v|² = 64 throughout:

```
+-1^64        2^28 sign vectors from a [64,28,16] code
+-2^16 0^48   30828 supports x 2048 signs
+-4^4  0^60   10416 supports x   16 signs   (10416 = C(64,3)/4, an SQS(64))
+-8^1  0^63      64 supports x    2 signs
```

A class needs |⟨u,u′⟩| ≤ 64/3, i.e. ≤ 21.

- **Free, and rigorous:** the 64 lines {±8eᵢ} are pairwise orthogonal, hence a class. That
  alone gives τ(64+k) ≥ 331 737 984 + 128 + τ(k), a gain of **+128** in dimensions 65–69.
  Small, but it would be a record in five dimensions and costs nothing.
- **Better:** the weight-4 family gives classes of 64 lines too — 16 pairwise-disjoint
  supports (a parallel class of the SQS(64)) times 4 sign patterns per support, the patterns
  pairwise at Hamming distance exactly 2 modulo complementation. Disjoint parallel classes
  give disjoint classes, so with λ(k) of them the gain becomes 128·λ(k): **+2560** at
  dimension 69.
- **Speculative:** in the full-support family two vectors have inner product 64 − 2d, so a
  class is a set of codewords with all pairwise distances in [22, 42]. The Kerdock code K(6)
  is a (64, 2¹², 28) code with weights in {0, 28, 32, 36, 64}, all inside that window — 2048
  lines, gain +4096 per group — but it would have to sit inside the [64,28,16] code ERS
  actually use, which is **unverified**. A Caro–Wei estimate on a random-like code of size
  2²⁸ also lands near 2¹¹ lines, so ~2000 is the right order.

**None of this is claimed here**, because none of it competes with dimensions 70 and 71,
where the Γ₇₂ cross-sections give factors of 3.8 and 7.9. It is recorded so the avenue is
not lost: anyone wanting dimensions 65–68 should start from the first bullet, which is
rigorous and takes an afternoon.

## Files

```
README.md   this file
```

There is no separate script: every number above comes from
[`common/kpoint_lp.py`](../../../common/kpoint_lp.py) with the Γ₇₂ one-point marginals of
[`common/theta.py`](../../../common/theta.py), the same two modules that
[`../../improved/dim70-71-gamma72-cross-sections/derive.py`](../../improved/dim70-71-gamma72-cross-sections/derive.py)
uses. To reproduce the table, call `run3(72, 6218175600, 8, [1,2,3,4], M1G, G)` for each Gram
above and `run3(..., target=cell)` for the realizability witness.

## References

See [`../../../common/published.py`](../../../common/published.py).
