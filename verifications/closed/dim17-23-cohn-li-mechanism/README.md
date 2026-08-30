# CLOSED — the Cohn–Li mechanism, dimensions 17 through 23

**Nothing is claimed here.** This package records that a whole family of constructions has
been pushed to its exact ceiling, so nobody spends time on it again. The one dimension that
was open when this project started — 19 — was closed during it, at somebody else's value.

```
python scripts/verify19.py     # the structural facts about the dimension-19 Cayley graph
python scripts/ceilings.py     # the Delsarte/theta ceiling for every n in 17..23
python scripts/alpha512.py     # alpha(G_0) = 160, proven optimal  (needs OR-Tools CP-SAT)
```

## The mechanism

H. Cohn and A. Li, *Improved kissing numbers in seventeen through twenty-one dimensions*,
[arXiv:2411.04916](https://arxiv.org/abs/2411.04916), hold dimensions 17, 18, 20 and 21 in
Cohn's table; B. S. Ho, [arXiv:2603.10425](https://arxiv.org/abs/2603.10425), holds 19.

For a binary constant-weight-8, minimum-distance-8 code *C* of block length *n* there is a
kissing configuration of size

> 4·C(n,2) + 128·|C|

— the 4·C(n,2) permutations of (±2,±2,0,…,0), plus 128·|C| vectors putting an **odd** number
of minus signs on each codeword. (Flipping the parity of the sign count is the whole idea:
the "even" version gives the Leech cross-sections.) Their Table 2.1, from repeatedly
shortening S(5,8,24):

| n | 17 | 18 | 19 | 20 | 21 | 22 | 23 | 24 |
|---|---|---|---|---|---|---|---|---|
| \|C\| | 30 | 46 | 78 | 130 | 210 | 330 | 506 | 759 |
| base | 4384 | 6500 | **10668** | **17400** | **27720** | 43164 | 65780 | 98256 |

Extra vectors v_c = (−1)^{c_i}·√(8/n) may be adjoined for any *c* orthogonal to *C* mod 2,
and two of them are compatible iff their Hamming distance is at least n/4, because
⟨v,v′⟩ = (n−2d)·8/n ≤ 4.

**So the extra term is exactly the size of a minimum-distance-⌈n/4⌉ code inside the dual of
C** — a punctured extended Golay code — which is a **maximum independent set in a Cayley
graph on F₂¹²** with connection set S = the non-zero words of punctured weight below n/4.
Once you see that, the ceiling is a few seconds of linear programming.

## The ceilings, and the verdicts

For an abelian Cayley graph the Delsarte LP bound equals the Lovász theta function, and both
are computed exactly here:

| n | d = ⌈n/4⌉ | \|S\| | components | θ per component | ceiling on the extra code | needed to beat the record | verdict |
|---|---|---|---|---|---|---|---|
| 18 | 5 | 51 | 1 × 4096 | **768.0000** | 768 | 1155 | **DEAD, provably** |
| 19 | 5 | 21 | 4 × 1024 | 384.0000 | 1536 | 1281 | **DEAD** — see below, α = 320 exactly |
| 20 | 5 | 5 | 128 × 32 | **16** | **2048** | 2049 | **DEAD, provably** |
| 21 | 6 | 21 | 4 × 1024 | **512.0000** | **2048** | 2049 | **DEAD, provably** |
| 22, 23 | | | | | ≤ 4096 total | gaps of 6733 and 27371 | **DEAD** |

**Dimensions 20 and 21 are exactly optimal for their own mechanism.** Dimension 20 is the
prettiest: |S| = 5 and those five words are linearly independent, so each component is the
5-cube Q₅, which is bipartite with independence number exactly 16 — hence 128 × 16 = 2048,
*precisely* Cohn–Li's gain. Dimension 21 gives θ = 512.0000 per component, again exactly
their 2048. Neither can be improved by this route at all.

**Dimension 18** needs a length-18 distance-5 code of size ≥ 1155 inside the dual. The
ambient is the 6-punctured Golay (4096 words) with |S| = 51, one component; the
ambient-restricted Delsarte/theta bound is exactly **768.0000**, and greedy plus (1,2)-swap
local search over 150 restarts **attains 768**. So α = 768 exactly, and 6500 + 768 = 7268
against the record 7654 — a factor 1.5 short, unbridgeable. (The unrestricted bound
A(18,5) ≤ 1289.5 was far too weak to decide this; the ambient restriction is what kills it.)

**Dimension 17** uses a different construction (the odd 16-dimensional configuration) and is
not comparable — the compatibility of v with the odd-sign vectors needs n ≥ 18. It is
handled separately in [`../dim09-19-record-maximality/`](../dim09-19-record-maximality/).

## Dimension 19: closed at exactly Ho's value

This was the one live dimension in the whole range, with α ∈ [320, 348] — room for up to
+112, i.e. τ(19) up to 12 060. It is now closed at **α = 320 exactly**, so
**τ(19) = 10668 + 4·320 = 11948 is the exact ceiling of this construction** and Ho's value
attains it. The chain, each link either enumerated or solved to proven optimality:

1. τ(19) = 10668 + 4·α(Cay(K,S)); the Cayley graph on the punctured Golay dual has four
   components, each a translate of Cay(K,S) with K = F₂¹⁰, |S| = 21.
2. **S has exactly one odd-weight element** *t* (the weight-3 word) and twenty even ones, so
   K = K₀ ∪ (K₀ + t) with K₀ the even-weight hyperplane. Both halves induce
   G₀ = Cay(K₀, S₀) (512 vertices, 20-regular) and the cross edges are a **perfect matching**
   x ↔ x + t, because x + y lies in K₀ and *t* is the only element of S outside it. Verified
   vertex by vertex by `scripts/verify19.py`.
3. Hence α(Cay(K,S)) = max{|A| + |B| : A, B independent in G₀, disjoint} — the maximum
   induced bipartite subgraph of G₀ — which is at most 2·α(G₀).
4. **α(G₀) = 160, proven OPTIMAL** (`scripts/alpha512.py`, CP-SAT). The model is exact:
   x₀ = 1 is free by vertex-transitivity, maximality holds for any *maximum* independent set,
   and each of the 40 000 girth-5 cuts (sum over a 5-cycle ≤ 2) is valid because an
   independent set in C₅ has at most two vertices. **Cross-checked** by re-solving with the
   maximality constraint dropped (`scripts/check_nomax.log`): still 160, still OPTIMAL, so
   the closure does not rest on that link. Without the girth-5 cuts the solver stalls at an
   upper bound of 348 — the cuts, not more time, were what mattered.
5. Therefore α(Cay(K,S)) ≤ 320, and Ho attains 320. **α = 320 exactly.**

Two structural facts made it reachable: the five (−12)-eigenvalue characters give a Clebsch
quotient — which is exactly where Ho's 64 × 5 = 320 comes from — and all five 4-element
blocks share the sum z = 499, making G₀ a ⟨z⟩-double cover of a 256-vertex graph with
α = 80 and maxbip = 160, so the lift family reproduces 320 and provably cannot beat it.

### What did **not** work as an upper bound

Recorded so it is not retried:

- the **clique–coclique bound** is useless: the graph is triangle-free, ω = 2, giving only 512;
- the **Hoffman ratio bound** gives 391 on Cay(K,S) and 192 on G₀;
- the **Clebsch-class sumset relaxation** (`scripts/classlp.py`, with fᵢ(m) = min |X + Ĝᵢ|
  computed exactly for every m) returns only 256, because it permits every class to sit at
  16 — a profile the real problem forbids, since a class at 16 would have to be Vᵢ-periodic
  for all five labels at once, hence empty or all of N. Capturing that needs the
  periodicity, not just the cardinality.

## History of the lower bound in dimension 19

| when | value | who |
|---|---|---|
| 1971 | 10 668 | the base configuration alone |
| Nov 2024 | 11 692 | Cohn–Li [CL24], arXiv:2411.04916 |
| 2026 | **11 948** | Ho [Ho26], arXiv:2603.10425 — a nonlinear subcode found via a 5-coclique in the Clebsch graph |
| this work | — | **no improvement, and none is possible**: 11 948 is proven to be the exact ceiling of the construction |

Dimension 19's record landed while this project was in progress. That is worth stating
plainly as a caution: **Cohn's table moves.** Two entries (12 and 19) changed during this
work. Re-fetch before claiming anything.

## Files

```
README.md                 this file
scripts/verify19.py       the structural facts: 4 components, the perfect matching, G_0
scripts/ceilings.py       Delsarte/theta ceilings for n = 17..23
scripts/alpha512.py       alpha(G_0) = 160, proven optimal by CP-SAT with girth-5 cuts
scripts/alpha512.log      its output
scripts/check_nomax.log   the same, with the maximality constraint dropped: still 160
scripts/lpbound.py        the two-point LP used for the ceilings
scripts/classlp.py        the Clebsch-class sumset relaxation that returns only 256
scripts/structure.py      the eigenvalue/character decomposition of G_0
scripts/ho.py             Ho's 320-vertex construction, reproduced
```

`alpha512.py` needs OR-Tools (`pip install ortools`); the others need only numpy and scipy.

## References

See [`../../../common/published.py`](../../../common/published.py).
