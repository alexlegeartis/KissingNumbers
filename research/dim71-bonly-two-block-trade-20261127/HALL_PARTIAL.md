# Addendum: Hall obstruction for partial 71D B-bank trades

Use the fixed coordinate model, full blocker definitions and exact pinned old 71D configuration in [the preceding proof](./README.md). This extends the complete-bank bound to **arbitrary subsets of specified banks**. It also proves independently that the small norm-eight first cut is not required to block append-only insertions in this class.

Define \(Y_b=\{(0,u,v):u,v\in F_b\}\) for the 4,095 nonzero B-only residues. The canonical sign action \(G=(C_2)^2\) flips the signs of the last two Leech blocks independently. It acts freely within each \(Y_b\) and on every fixed-first-block blocker \(z=(\lambda,z_2,z_3)\) below, since neither tail block vanishes. Collapse each set to orbits of size four. Distinct old orbits have disjoint sets of actual old points.

## Finite certificate

The 6+6 displacement family gives a blocker of an orbit on each ordered axis combination for which the two local absolute products add to 6 or 8. For a fixed candidate orientation, choice of the old sign in each coordinate is unique within any resulting old orbit: for \(w_i=z_i-y_i\) of norm six and \(y_i\) of norm eight,
\[
\langle z_i,y_i\rangle=(\|z_i\|^2+2)/2>0.
\]
Thus raw degree per fixed candidate and orbit degree coincide. Further, in any fixed bank label \(b\), each old 6+6 blocker orbit meets at most one candidate orbit. This follows from the fixed signed old-coordinate uniqueness in the main result, then equivariance under the same sign action.

With all 4,095 banks included, the exact replay gives:

* For **every odd bank**, all 576 unsigned ordered pairs have at least 13 old 6+6 blocker orbits.
* For **even banks**, every pair with any 6+6 blocker has at least 4. There are precisely 1,248 unsigned ordered pairs with no 6+6 blocker: 240 banks have four each and 18 have 16 each. The remaining 1,789 even banks have no such pair.
* Every one of these 1,248 even residual pairs has at least 40 distinct old 4+8 or 8+4 blocker orbits. Each of those old orbits meets at most two residual pairs within any fixed even bank. This counts distinct literal pairs of canonical signed-block coordinate vectors, not displacement presentations, so these are true neighbor degrees.

The right signatures that separate the two pools are also elementary. Even 6+6 blockers have tail norms \((10,6),(6,10)\) or \((6,6)\). Odd 6+6 blockers and even 4+8 blockers have no norm-six tail coordinate: their possible signatures are \((12,4),(8,8),(4,12),(8,4),(4,8)\). These follow from the tail norm equations \(8+6-2h\) and \(8+4-2h,8+8-2j\), with \(h+j=6\) or 8 and the parity/upper bounds on local products. All claimed witnesses are part of the actual retained old stock, as in the preceding proof.

## Hall argument

Choose any \(m_e\) even and \(m_o\) odd bank labels subject to
\[
m_e\le4,\qquad \frac{m_e}{20}+\frac{m_o}{13}\le1.
\]
Include all candidate orbits in these banks initially. Each sends one unit uniformly across the indicated blocker orbits: use 6+6 if available, 4+8/8+4 for the even residues with none. Old blocker orbits with a norm-six tail coordinate receive at most \(m_e/4\), since at most one orbit per even bank can charge them. Other old blocker orbits receive at most \(2m_e/40+m_o/13\), since at most two residuals per even bank or one per odd bank can charge them. Thus each old orbit receives at most one unit.

Any collection \(S\) of candidate orbits sends \(|S|\) charge to its old neighbors, proving the finite Hall condition. There is therefore an injective choice of old blocker orbits. For each matched edge choose an actual conflicting pair \((y,z)\) of representatives, and map \(gy\mapsto gz\) for \(g\in G\). Since both actions are free, this lifts to an injection from individual candidate points to individually conflicting old points. It applies to any subsequently selected partial subset of the candidate points, without requiring sign balance or whole groups. Arbitrary deletions of shared old blockers cannot turn such a trade positive.

In particular this excludes gain from arbitrary partial trades on any set of up to **four** banks, as well as e.g. 13 odd banks or 4 even plus 10 odd. This does not cover arbitrary partial trades with larger support outside the stated inequality. No bound changed.

The [standalone companion](verify_partial.mjs) uses the same three hash-pinned inputs as the main checker and reconstructs the exact orbit degrees in both pools. Run `node research/dim71-bonly-two-block-trade-20261127/verify_partial.mjs .` from the repository root.
