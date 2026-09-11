# `verifications/` — four tiers

Every construction this project produced lives in exactly one of four directories,
according to what is currently true about it. The distinction is the point of the layout:
a repository that only shows its wins is much less useful to the next person than one that
shows where the doors are shut.

## [`improved/`](improved/) — 47 dimensions, currently the best known

Nine packages, one per **construction idea** rather than one per dimension, because several
ideas each cover a range. Each has a README with the idea in prose, a description of every
file, and the history of the lower bound in its dimensions.

| package | dims | status |
|---|---|---|
| [`dim25-lens-heads`](improved/dim25-lens-heads/) | 25 | apparently new — joint work in progress with H. Cohn and B. Lindow; 1006 lens heads, one non-lattice equator point |
| [`dim26-27-iota-triangles`](improved/dim26-27-iota-triangles/) | 26, 27 | apparently new — joint work in progress with H. Cohn and B. Lindow; the coset triangle on every triangle of directions, the side chosen per head |
| [`dim38-leech-large-codimension`](improved/dim38-leech-large-codimension/) | 38 | apparently new |
| [`dim39-ers-constant-weight`](improved/dim39-ers-constant-weight/) | 39 | apparently new — but check Sun–Wang's latest revision |
| [`dim49-63-p48-caps`](improved/dim49-63-p48-caps/) | 49–63 | apparently new — 62 and 63 are reached here too, but claimed from the chain below |
| [`dim62-63-ers-chain`](improved/dim62-63-ers-chain/) | 62, 63 | apparently new — the Edel–Rains–Sloane chain, which beats the caps in both |
| [`dim68-69-gamma72-cross-sections`](improved/dim68-69-gamma72-cross-sections/) | 68, 69 | apparently new — the k-point LP with the Gram exhibited |
| [`dim70-71-gamma72-cross-sections`](improved/dim70-71-gamma72-cross-sections/) | 70, 71 | apparently new — **the claim most in need of a literature check** |
| [`dim73-95-gamma72-caps`](improved/dim73-95-gamma72-caps/) | 73–95 | apparently new |

## [`recovered/`](recovered/) — already published, re-derived here

One package, and **nothing in it is claimed**. Both of its values were in the literature
before this project reached them, and both were reached here independently, by a method that
knows nothing of either. That is exactly why the package is kept: reproducing a published
count from an independent route is the sharpest check the cross-section machinery has.

| package | dims | who has it |
|---|---|---|
| [`dim46-47-p48-cross-sections`](recovered/dim46-47-p48-cross-sections/) | 46, 47 | **47**: Boyvalenkov–Cherkashin, *Results in Mathematics* **80** (2025), Paper No. 3; equation (3) of [arXiv:2312.05121](https://arxiv.org/abs/2312.05121). **46**: Sun–Wang, [arXiv:2607.20359](https://arxiv.org/abs/2607.20359), Table 3; the ingredients were in Ozeki's 2016 Siegel theta tables |

Cohn's table has absorbed neither, still recording 9 741 412 and 5 318 060. That is a fact
about the table, not a result of this project, and it is not counted as one.

## [`superseded/`](superseded/) — was a record, is not any more

Three packages. None is wrong; all were simply beaten, one by other people and two by this
project's own later work.

| package | dims | beaten by |
|---|---|---|
| [`dim44-45-p48-cross-sections`](superseded/dim44-45-p48-cross-sections/) | 44, 45 | Sun–Wang, arXiv:2607.20359v3 — **and their values land inside this project's LP brackets**, which is the sharpest available cross-check on machinery four surviving claims depend on |
| [`dim38-ers-constant-weight`](superseded/dim38-ers-constant-weight/) | 38 | this project's own Leech cap construction, by +21 376 |
| [`dim25-cap-level`](superseded/dim25-cap-level/) | 25 | this project's own lens-head configuration, by +511 |

## [`closed/`](closed/) — mechanisms at their exact ceiling

Ten packages. Nine carry no claim: they are the negative results, kept because they are what
stops the next attempt from wasting time, and several are *proofs of optimality for a
mechanism* rather than failed searches. The tenth is the exception and is marked as such
below — it settles dimension 96 **against** the cap construction, and the claim it leaves
standing is Edel–Rains–Sloane's rather than this project's.

| package | what it settles |
|---|---|
| [`dim17-23-cohn-li-mechanism`](closed/dim17-23-cohn-li-mechanism/) | the Cohn–Li construction reduces to a Cayley-graph independence number; **provably exhausted in 18, 20, 21**, and **dimension 19 closed at exactly Ho's 11 948** |
| [`dim22-23-maximal-cross-sections`](closed/dim22-23-maximal-cross-sections/) | Λ₂₁, Λ₂₂, Λ₂₃ are **maximal spherical codes** with exact min-max cosines √(8/29), √(3/11), √(4/15); Leech level sets cap at exactly 93 150 by a Hall/König argument |
| [`dim09-19-record-maximality`](closed/dim09-19-record-maximality/) | **every published record in dimensions 9–19 is maximal**, with exact margins; and layering over a record equator is measured and shown to be always a mistake |
| [`dim32-44-ers-audit`](closed/dim32-44-ers-audit/) | dimensions 32–44 are **exactly at** the Edel–Rains–Sloane value; **Cohn's table is stale in 32, 33, 34, 37**; and exactly where a code improvement would pay |
| [`dim37-cap-shortfall`](closed/dim37-cap-shortfall/) | the cap construction at k = 13 falls 4% short, measured end to end; it needs a better *class*, not a better cover |
| [`dim69-gamma72-three-point`](closed/dim69-gamma72-three-point/) | dimension 69 does not reach its floor; **one large tempting number discarded for want of a realizability proof**; also dims 53–55 and a recorded partial win in 65–68 |
| [`class-problem-upper-bound`](closed/class-problem-upper-bound/) | the central open computational question: 248 lines against a 425.45 bound that survives every two-point method, the three-point SDP and the subconstituent split |
| [`antipode-construction-ceiling`](closed/antipode-construction-ceiling/) | the antipode construction **is** a max-weight clique problem, whose arms lie in a 30° cap — a structural reason it cannot give low-dimensional records |
| [`dim96-ers-takeover`](closed/dim96-ers-takeover/) | **the one claim in this tier.** At 96 the chain (96, 24, 6, 1) reaches 12 886 999 232 against the 6 480 558 568 the layered construction manages, so the cap route loses outright; also the exposure measurement for 93–95 |
| [`dim17-layered-family`](closed/dim17-layered-family/) | τ(17) = 5346 + 2α with **α = 192 exactly**, so τ(17) = 5730 and there is no dimension-17 record in this family. This was the repository's one live lead until `K'` was noticed to be closed under XOR: `Cay(K', W_4)` is then a translation scheme, Delsarte's LP gives 192 on the nose, and its dual is six rational numbers verified in integer arithmetic — replacing a CP-SAT bound that ranged over 192, 193, 195 and 238 across runs |


## Reading order, if you are picking this up

1. [`../README.md`](../README.md) and [`../RESULTS.md`](../RESULTS.md) — what is claimed.
2. `python ../audit.py` — that the claims are mutually consistent.
3. [`superseded/dim27-triple-partition/`](superseded/dim27-triple-partition/) — the smallest
   complete package, verifiable from a coordinate file alone with no third-party packages.
   It is the model the others aim at, and since 2026-09-08 it is superseded in its own
   dimension by [`improved/dim26-27-iota-triangles/`](improved/dim26-27-iota-triangles/).
4. [`../common/PROOF-kpoint.md`](../common/PROOF-kpoint.md) — the mathematics behind the four
   cross-section claims, dimensions 68 to 71, and behind the recovered 46 and 47.
5. `closed/` — before starting anything, to find out whether it is already known to be dead.
6. [`closed/dim17-layered-family/`](closed/dim17-layered-family/) — the best worked example of
   the repository's most useful habit: when a search will not converge, stop searching and ask
   what structure the object has. The answer was "it is a group".
7. [`../KNOWLEDGE.md`](../KNOWLEDGE.md) — 127 sections of working notes, most of
   them about things that did not work.
