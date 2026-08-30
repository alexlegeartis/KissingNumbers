# CLOSED — the antipode construction is a max-weight clique problem, and that caps it

**Nothing is claimed here.** This records a structural theorem: the antipode construction —
the method behind Sun–Wang's dimensions 39, 43 and 45 — reduces exactly to a max-weight
clique problem, and the reduction proves it **cannot produce low-dimensional records**. This
is an obstruction, not a search failure, so the avenue does not need retrying.

```
python antipode3.py         # the measurement over the Leech, k = 1..12
```

## The reduction

Sun–Wang ([arXiv:2607.20359v3](https://arxiv.org/abs/2607.20359)) and Chen et al. (2025)
state the construction as: for a lattice *L* of minimal norm μ, a rank-*k* cross-section with
Gram *K*, and *M* = *K*⁻¹, pick a finite S ⊂ *M* whose maximal pairwise squared distance δ is
< μ. Then

```
N_i = sum over j with |u_i - u_j|^2 = delta  of  fibre(u_j - u_i),
fibre(c) = #{ w in L : |w|^2 = mu, U(w) = c },   U(w) = (<u_1,w>, ..., <u_k,w>)
```

and the packing's kissing number is max_i N_i. Chen et al.'s Theorem 2 takes S to be the
dual-basis simplex {0, m₁, …, m_k} — one point of a large search space. The whole space
collapses:

> **Theorem.** Over all admissible S,
>
> ```
> max_S max_i N_i  =  max over delta in (0, mu) of MaxWeightClique(V_delta, E_delta, fibre)
> ```
>
> where V_δ = {v ∈ M : |v|² = δ, fibre(v) > 0} and v ~ v′ iff ⟨v,v′⟩ ≥ δ/2.

*Proof.* (≤) Fix S and the index *i* attaining the max; let v_j = u_j − u_i run over the
diameter partners of u_i. Each |v_j|² = δ, and for two partners
|v_j − v_l|² = |u_j − u_l|² ≤ δ; with |v_j|² = |v_l|² = δ that is exactly ⟨v_j, v_l⟩ ≥ δ/2.
So {v_j} is a clique of weight N_i. (≥) Given a clique V, set S = {0} ∪ V. Then |v|² = δ and
|v − v′|² = 2δ − 2⟨v,v′⟩ ≤ δ, so the maximum pairwise distance is exactly δ, every pair (0,v)
is a diameter pair, and N₀ = Σ_{v ∈ V} fibre(v). Centring at 0 is free because *M* is a
lattice and only the difference set matters. ∎

## Corollary: the cap ceiling

⟨v,v′⟩ ≥ δ/2 with |v| = |v′| = √δ says the angle between any two arms is at most 60°, so
**all arms lie in a cap of angular radius 60° about any one of them**. So the antipode
kissing number is at most the fibre mass of the heaviest such cap of a single sphere
|v|² = δ — it can never reach even half the shell.

**An earlier version of this section said 30°, and that is wrong.** Pairwise 60° does not
give a 30° cap: two arms lie within 30° of their bisector, but three arms pairwise at
exactly 60° have centroid cosine 2/√6 and need arccos√(2/3) = 35.26°. Nothing measured
below depended on it — `antipode3.py` computes max-weight cliques directly and never uses a
cap radius — and the corollary that actually bounds the star, that the clique is equilateral
of rank at most k, is the next one. In particular the antipode beats the bare cross-section fibre(0) only when
the fibre distribution is very flat.

**Corollary (why the simplex is the worst case).** If δ is the minimal norm of *M* then
v − v′ is a non-zero *M*-vector of norm ≤ δ, hence *= δ*: the clique is equilateral, its Gram
is (δ/2)(I + J), positive definite of rank *m*, so *m* ≤ *k*. That is precisely the
dual-basis simplex regime. Larger δ allows larger stars but lighter fibres.

## The measurement

`antipode3.py`, over the Leech lattice with greedy cross-sections, k = 1…12. The optimum is
always attained at small δ with a 1–4 arm star, and always loses:

| dim | bare cross-section | best star | Λ_dim | record |
|---|---|---|---|---|
| 23 | 93 150 | 47 104 (1 arm, δ = 1/4) | 93 150 | 93 150 |
| 22 | 49 896 | 41 472 (2 arms, δ = 1/3) | 49 896 | 49 896 |
| 20 | 15 540 | 7 040 (2 arms, δ = 3/5) | 17 400 | 19 448 |
| 18 | 5 820 | 4 440 (4 arms, δ = 11/15) | 7 398 | 7 654 |
| 16 | 2 640 | 2 100 (4 arms, δ = 4/5) | 4 320 | 4 320 |
| 12 | 624 | 396 (3 arms, δ = 1) | 648 | 841 |

**So the antipode is a ~10% lever that only pays when the bare cross-section is already
within about 10% of the target** — which is the dimension 44-to-47 situation, where Sun–Wang
use it successfully, and *not* the low-dimensional one, where the bare cross-section is
8–36% below the record.

## Two side facts recorded while doing this

- **The chain u₁ … u_k must NOT be built from mutually orthogonal minimal vectors.** Only the
  cross-section has to be orthogonal to all of them; restricting the candidates to the
  surviving orthogonal set gives 43 164 in dimension 22 instead of Λ₂₂'s 49 896. Greedy over
  the whole shell recovers the laminated values 93150 / 49896 / 27720 / 17400 / 10668 / 7398
  exactly (the chain Z, A₂, A₃, D₄, D₅, E₆, E₇, E₈).
- **The antipode packing is |S| translates of the cross-section lattice** L′ = L ∩ K⊥ at
  minimum distance √(μ − δ): P = {π(w) : w ∈ L, U(w) ∈ S}, and
  |π(w) − π(w′)|² = |w − w′|² − |U(w) − U(w′)|² ≥ μ − δ. Contacts *inside* one translate are
  at √μ > √(μ − δ) and therefore never count. Getting them to count would need the translates
  to be √μ apart mod L′, i.e. covering radius ≥ √μ; the Leech's is √2, so that door is shut.

## Files

```
README.md              this file
antipode3.py           the reduction and the k = 1..12 measurement over the Leech
antipode_leech.py      the fibre computation
```

## References

See [`../../../common/published.py`](../../../common/published.py). Sun–Wang's own successful
use of this construction, in dimensions 44 and 45, is discussed in
[`../../superseded/dim44-45-p48-cross-sections/`](../../superseded/dim44-45-p48-cross-sections/).
