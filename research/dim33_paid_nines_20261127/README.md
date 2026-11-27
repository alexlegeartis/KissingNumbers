# Lower90: actual 33D paid octad trades and anisotropic nine-site layer

**Scoped research note; no new lower bound.** Repo's independently replayable count remains **370,528** ([current witness and verifier](https://github.com/alexlegeartis/KissingNumbers/tree/main/verifications/improved/dim33-joint-signed-refill)), composed of 131,072 flat binary signs, 1,837 specified octad banks of 128, and 4,320 explicit rows (2,216 sign octads, 172 stars, 1,932 roots). A separate dataset by P. Carder reports a greater **370,568**, [dataset DOI](https://doi.org/10.5281/zenodo.24577265); we have not independently replayed Carder's coordinates. Qihong Lin's [published 370,240](https://arxiv.org/abs/2611.03682) is no longer the actionable target. Strict target is 370,569.

## Actual price of a new complete bank exchange

The candidates were cross-checked directly against all 1,837 full carrier masks and every one of the current 4,320 explicit rows, not the older 370,408 support source alone. The nine supplied triples each remove one full carrier and add two whose mutual intersection is at most four; each new carrier touches no other old full carrier in five or more positions. Counts below report the number of *distinct* explicit rows that must be removed for choosing the stated complete 128-sign parity family on each entrant. They are identical for all four parity assignments.

| Old mask | New mask A | New mask B | Forced explicit deletions |
| ---: | ---: | ---: | ---: |
| 1074431249 | 285933841 | 1074955537 | 336 |
| 4580344841 | 4580183809 | 4580245519 | 612 |
| 4295055681 | 4312861953 | 4563472593 | 744 |
| 4295055681 | 4564520193 | 4563472593 | 768 |
| 4429259523 | 4514185473 | 4563470097 | 584 |
| 4429259523 | 4546691329 | 4563470097 | 584 |
| 4568654401 | 4581257217 | 4564453137 | 516 |
| 2170621961 | 4578152449 | 285282571 | 664 |
| 2170621961 | 4578152449 | 2969637121 | 600 |

The nominal gain 128 is erased even by the best of these trades (result 370,320). This is an exact conclusion for these tested complete sign trades, not all support or signed trades. For future block releases, a clean exact accounting mechanism is signed 5-face capacity: with five common coordinates fixed at ±2, two norm-32 rays have tail dot ≤−4 and each tail has squared norm 12; the square of the sum gives ≤4 rows per signed face. A full parity octad already fills all four slots on each of its 32 signed faces. A fixed unsigned 5-face therefore caps its own complete-octad replacements at 128 even if signs are arbitrarily reassigned. The 1→2 trades above use different overlapping faces of the released carrier, so this individual charge does not exclude them; the actual sparse prices do.

## A structurally different nine-support move, and its exact prices

The complete monotone enumeration of uncovered coordinate sets (no five-set contained in a fixed complete old octad) has cardinalities 134,464, 79,706, 3,044, 197, 18, 0 for sizes 5 through 10. In particular all **18** uncovered 9-sets contain the unique zero coordinate of the dense bank, physical index 32; there are no larger uncovered sets. This explains a genuinely new layered direction: start with 9-site rays, putting a signed a=√(32/7) at coordinate 32 and signed b=√(24/7) at each of the other eight sites. The exact squared norm is 32. Every ray in this 9,216-element reservoir obeys all dense bounds since 8b<16, all root bounds since a+b<4 and 2b<4, and all complete octad bounds since the support intersections have size≤4 and a+3b<8. The corresponding uniform a=b=√32/3 carrier also meets those stocks.

Between two anisotropic rays, let q be signed overlap on their ordinary active coordinates (integer). For equal polar signs the bound is q≤3; for opposite polar signs it is **q≤6**, one unit looser than the uniform 9-site radius profile's q≤5. These equalities are immediate from raw dot ±32/7+24q/7. Every signed pair with an old explicit row can be checked wholly with integer linear combinations of rigorous rational enclosures of a,b. Exactly the same explicit rows obstruct either the anisotropic or the uniform profile, **8 to 32 old rows per individual new candidate**, all literal old octads. The union has 1,864 old signed points; none of the 9,216 candidates is unpaid. On any single support, each polar sign class gives a length-eight binary code of distance≥3, with at most 20 words (append the parity bit and apply [Best–Brouwer–MacWilliams–Odlyzko–Sloane 1978, A(9,4)=20](https://experts.umn.edu/en/publications/bounds-for-binary-codes-of-length-less-than-25/)). Thus a lone carrier would already be too small even before paying or considering the opposite-polar interface; a gain needs coupled carriers.

The same price obstruction applies to **any strictly positive amplitudes≤2 on any one of these 18 support/sign patterns**: write w_j=2(1−d_j), with nine entries, 0≤d_j<1 and Σ_j w_j²=32. Then Σ_j(2d_j−d_j²)=1 and Σ_j d_j<1. Any explicit octad blocker of the uniform profile has signed overlap q≥5. Its new raw dot is 4(q−Σ_j ε_j d_j)>16 because ε_j∈{−1,0,1}. All its original old blocker rows are still required payments irrespective of the new amplitude profile. This lemma does not cover amplitudes>2 (including the anisotropic profile, checked separately) or other supports. It is not a many-donor ceiling, because a blocker may pay for several simultaneous donors.

The 18 supports assemble into six natural Hamming-overlap components with ≥8 intersections, each on only ten coordinates. Their carrier counts are 3,4,3,2,4,2. This makes a future mathematical exchange test concrete: solve signed capacity jointly in these ten-coordinate components against their actual intersections with the 27 affected old octad carriers, using the fixed-face slack when nine-site rows replace eight-site rows. Half of the logic is not an arbitrary global ReduMIS cutoff: same-pole conflict is entirely a code condition (q≤3), cross-pole is q≤6, and prices are fixed signed 5-/6-/7-/8-face incidences. No jointly improving witnessed code was found, and no optimality claim for this portfolio.

Replay: [`verify.mjs`](verify.mjs) needs the existing 33D joint-signed-refill witness.npz (SHA256 90d130ba1bb96172fae7b7ae10d1169747b76f2e720946bc2350c99a08d427c9). Node ≥20: node verify.mjs witness.npz. It hashes and decodes literal masks, verifies all covered/uncovered sets and all prices, and uses no solver or floating point angle checks.
