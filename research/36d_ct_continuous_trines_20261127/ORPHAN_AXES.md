# 36D held-out heavy bodies: exact exclusion of the natural mixed ternary axes

Same pinned live 506,628-point CT lower-bound witness and Yang et al. outside 484,760 comparison as the Lower79 radial proof. The current [neutral-reheat ledger](https://github.com/alexlegeartis/KissingNumbers/blob/a304d9f04756c67f07c7df9fc846b989c1bb0826/verifications/improved/dim36-ct-neutral-reheat/witness.json) removes all three rows of 1,950 heavy trines. It restores all three rows on 816 heavy bodies, leaving **1,134 genuinely unoccupied heavy rays** \(v\) of norm squared eight (one fixed physical block). The live triple-root sector comprises 106,326 rows. The exact replay reconstructs and checks both counts through the actual insertions/deletions and all three ancestor files, not the unswitched parent.

For a held-out heavy ray \(v\), take a phase \(r\) from the old straight root trine in its companion block and \(s\) from the root trine of the opposite label in its spare block. Both have raw norm squared four. Consider every real mixed row
\[
  (v,a r,b s),\qquad a^2+b^2=1,\quad ab\ne0,
\]
over all nine phase pairs, so all four sign choices of the magnitudes are permitted. This is an exact no-append result in that three-block family, not a classification of arbitrary angles in the two root planes.

The old \(6+6\) rows at body product \(5\) give tail scalar \(1\) or \(-2\) as extrema (respectively) against each of the two positively oriented root directions, for both companion blocks and all 1,134 heads. These extrema are checked and persist under phase rotation. Thus a negative coefficient must have magnitude at most \(1/2\). Both coefficients cannot be negative, since their squares must sum to one.

The exact live \(4+4+4\) sector supplies, for each head and each phase pair, old rows having head product 4 and projected tail normals in **all four** of these four desired cones:
- first coordinate at least 2 and second at least 1;
- first at least 1 and second at least 2;
- first at least 2 and second at most −1;
- first at most −1 and second at least 2.

Thus if \(a,b>0\), the first two old rows imply \(2a+b\le2\) and \(a+2b\le2\). But if \(a\ge b\), \(2a+b>2\) on the positive unit quadrant: as a function of \(a\in[1/\sqrt{2},1]\), \(2a+\sqrt{1-a^2}\) is strictly concave, with endpoints \(3/\sqrt2>2\) and 2. If \(b\ge a\), the other wall gives the same contradiction. If one coefficient is negative, the other positive coefficient has larger magnitude because the negative one is at most \(1/2\). Use the third or fourth wall accordingly to contradict \(2\max(|a|,|b|)+\min(|a|,|b|)>2\). This is not a check at fixed angles; it covers every real ratio on these signed axes.

The companion standalone Node checker `verify_orphan_axes.mjs` takes, in this order, pinned C = `dim36-ct-nonlocal/certificate.json`, T = `dim36-ct-triangle-refill/certificate.json`, D = `dim36-ct-double-trine/witness.json`, W = `dim36-ct-neutral-reheat/witness.json` under `verifications/improved/`. It recomputes both the six-layer pressure and the four triple-root wall flags on all 10,206 phase-paired rays. Current guarantee is solely against the full old code, with no new deletions. Genuine bypasses are rotating inside the two root planes rather than taking their six axial phases, moving the norm-eight head itself, or a coupled trade of multiple old triple-root blockers. No new record. No native or Git.
