# 21D: at least seventeen original holes must change in every improved line-tag construction

Lower87, 27 November 2026. No new kissing lower bound. The benchmark is Kevin Li's [primary Hexagon paper, §3](https://hexagonmath.org/2610.00248v1) ([literal PDF](https://hexagonmath.org/pdf/2610.00248v1)), not arXiv. The bound is 31,560 = 27,720 odd sparse points + 3,840 Golay signs with projective-line weights. This result only concerns local modifications of Li's actual holes and the fixed line-tag architecture.

In one of the four isomorphic 1,024-word components, write G for the additive span of the 21 lines. It is rank ten. The old omitted set is H=ker η (64 words), with tags exactly as in Li Table 2. Tags 1,2,3 occur three times apiece and the other 12 nonzero tags occur once. The three repeated tags are the nonzero elements of a 2D subgroup of the tag space. Shift the component so that the old holes are this H. The 21 lines all have weight five, so the full component graph is bipartite; H contains 32 holes of each parity.

**Claim.** No new dominating hole set of smaller size than 64 can differ from H by removing at most **sixteen** original holes in this component. This applies to arbitrary replacement holes (possibly nonindependent), and even to those that do not admit any compatible labeling. The necessity of domination follows from the local line rule: every occupied vertex has a missing neighboring word along its chosen line. So the claim is stronger than an exclusion only for independent dominating sets.

## Mandatory vertices and exact incidence weights

The old hole h has itself and twelve private neighboring vertices h+P, where P ranges over the singleton-tag lines. No existing hole other than h covers any of these twelve leaves, and no two old holes are adjacent because no line has tag zero. If h is removed, new holes must cover all thirteen mandatory vertices. The sets of old private leaves are disjoint between their owners.

Write k_p for the count of removed old holes of parity p. Let s_p and t_p be counts of new holes of respectively singleton tag and triple tag at parity p. To cover the removed holes themselves, the new holes of opposite parity must satisfy

```
k_p ≤ s_(1-p) + 3 t_(1-p).
```

A new hole of opposite parity covers a private leaf of parity 1−p only if it is literally that leaf. Thus it contributes at most one such charge when its tag is singleton, zero otherwise. A new hole x of parity p may cover several private leaves from the removed holes of its own parity. Define its 32-entry charge vector on the parity-p old holes by

```
a_h(x) = number of pairs (P,Q) with P a singleton-tag line,
         Q any line, and x = h + P + Q.
```

Distinct unordered sums of two lines are distinct: in any proposed relation of four distinct projective lines, a point on one line outside the other three gives a contradiction. Thus a_h(x) ≤2; multiplicity two means both lines in the sum are singleton-tag. The exact per-type census is:
- singleton x: four old owners have charge 2 and nine charge 1;
- triple x: six old owners have charge 2 and none charge 1.

For the second statement, adding a repeated tag preserves the set of twelve singleton tags, so all twelve ordered pairs give six unordered pairs. For the first, adding a singleton tag maps eight singleton tags to singleton tags (four unordered pairs), three to repeated tags (nine ordered contributions), and one to zero (none). The literal checker verifies this rowwise. For any chosen k old owners of a fixed parity, one new singleton can therefore cover at most min(2k,k+4,17) leaves, and one triple at most min(2k,12). Private-leaf coverage requires at least 12k_p on side p, including the at most s_(1−p) from opposite parity.

## Coupled capacity certificate

For any two or three distinct same-parity new holes of prescribed types, add their 32 charge vectors. The largest possible sum on k removed old holes is the sum of the k largest entries. Tabulate the maximum of this quantity for SS,ST,TT and SSS,SST,STT,TTT. These maxima bound simultaneous charges even though some leaves might be counted repeatedly; repeated counts can only weaken the inequality. Exact computation at one parity suffices because translating by an odd h∈H maps between the parities. To compute all triple maxima without repeating symmetrical triples, normalize one new hole to the representative of its η-tag under the 32-element even subspace of H; the routine checks 1,717,215 triples, with no numerical calculations.

For candidate counts (s_p,t_p), the total charge can be bounded only from the largest three charges of each type: all remaining members of that type are at most its third largest. Enumerate these ordered integer sequences, checking every applicable pair/triple cap on the largest values. This is an upper relaxation, so exclusion remains rigorous. Enumerate 1≤j<k≤16, k0+k1=k and s0+t0+s1+t1=j. A j=0 exchange fails immediately because removed holes otherwise have no remaining hole neighbor. After requiring removed-owner domination and all 12k_p private charges, single-donor inequalities alone exclude through k=13, pair inequalities through k=15, and the mixed triple inequalities through **k=16**. Twelve count-level cases remain at k=17; they are not claimed physically realizable or impossible. See `lower87_21d_dominating_16_verify_20261127.mjs`, which recomputes every capacity and integer comparison from literal Table 2 lines.

For intuition, after only unary capacities the k=14 case is forced to have twelve removed holes on one parity, nine same-parity singleton new holes and four opposite-parity triple new holes, up to parity swap. All nine singleton donors would have to saturate sixteen private leaves on exactly the same twelve old holes. The independent direct check found all 3,456 such size-twelve support masks distinct (each is one donor's thirteen nonzero old owners minus one of its nine charge-one owners). The pair comparison also eliminates this case.

## Order-four odd quotient, a separate local scope

The previously known [Lower28 projective quotient audit] already excluded all even translation subgroups of order four or more and all Singer order-seven orbit-union covers of at most 63 holes, without assuming independence. Its durable notes are `/user_notes/long_form/kissing_lower_push28_20261125.md` and `/agent_notes/research_checkpoints/2026-11-25-kissing_lower28_projective_replay.md`. The recent Singer-order-seven independent-domination calculation is therefore a narrower alternate replay; this round does not claim it as new. The order-21 replay, if discussed, is a separate stated symmetry and not assessed here for novelty.

A subgroup W of order four containing odd words has even subgroup only of order two, so that older even-subgroup exclusion does not apply. Restrict here to the 496 such W≤H. Project to the 256-word quotient, with sixteen original hole cells H/W. The thirteen mandatory vertices around each removed old cell remain distinct: their η-tags are zero and the twelve different singleton labels, and W≤ker η. For each actual new quotient cell list exactly which mandatory sites it covers around each old one. No pair of new quotient cells can cover all thirteen around even one original cell. For three new cells, translate a removed old one to zero; since no pair covers it, each of the three covers a nonempty portion of its mandatory set. Enumerate triples using the exact 13-bit coverage masks; 222 triples cover the distinguished old cell across the 496 quotients. They cover at most one original cell's whole mandatory set in every case. Thus no 2→1, 3→2, or 4→3 replacement of old cells, even allowing nonindependence and ignoring labels. This does not exclude larger whole-cell trades, W not contained in H, or arbitrary period-two partial replacements. Replay `lower87_21d_odd_quotient_local_verify_20261127.mjs`.

The useful next positive route is not another handful of Li holes: the proof now forces at least seventeen changed holes in an improved component. Build a replacement on both parity sides simultaneously, or change the actual amplitude/carrier geometry. No native, repository, or external communication from worker.