# Midrange whole-layer campaign, 8 November 2026

**New verified lower bound: 251,200 in dimension 31**, four above the preceding
251,196-point repository construction. The portable package is
[`dim31-latitude-contact-exchange`](../../verifications/improved/dim31-latitude-contact-exchange/).
The frozen witness passes two independent complete exact pair enumerations,
including a direct-integer scan without antipodal compression.
All source files are the current complete 27–31D repository certificates;
source hashes appear in each report. No source or registered certificate is
modified. Search pools and numerical optimizer outputs are not geometric proofs.

## Exceptional tail-root replacement

`spinor_gate.py` tests a whole-axis replacement of the D6/D7 tail geometry
by the exceptional E6/E7 spinor roots. Every old point using the designated
tail axis is removed. Each distinct norm96 core body is paired with every
new spinor tail of squared norm32. All remaining old points have zero on
the radical coordinate, so the complete fixed-point gate uses exact integer
products, not a numerical radical tolerance.

- 30D, axis28: 1,485,440 candidates, none compatible with the retained code.
- 31D, axis24: 3,066,880 candidates, 256 compatible with the retained code,
  versus 16,260 mandatory removals. New/new compatibility was not needed to
  exclude improvement in this particular replacement.

A separate general limitation explains why changing the tail root system
alone is insufficient. For a fixed squared-norm96 core body, two norm32
 tails must have product at most -32. They are antipodes, so the body admits
at most two points. The current layers already nearly saturate that count.
This observation says nothing about adding new core bodies or changing radii.

## Complete tetrad-head grammar

`global_tetrads.py` generates all 10,626 core tetrads, every sign, and all
norm64 tail vectors of the forms signed four-coordinate tetrads and a single
coordinate ±8. Any such core body clears every old equatorial point, by the
checked bound four times the largest-four-coordinate l1 norm ≤64. All
remaining blockers are calculated explicitly. Every incumbent point is
removable; candidate selection does not freeze the other layers.

The four-blocker pool is a restriction, not a complete search over all
exchanges. A numerical MILP timeout is unresolved, not a certificate of no gain.
For 29D the complete 15,301,440-point grammar has 6,461 new candidates at four
blockers. The initial 60-second optimizer has no improving incumbent and a
numerical gain bound of 713, hence remains very far from closed.

### Exact capacity closure in 27D and 28D

`fiber_capacity.py` directly verifies six occupied tail directions with
48 points each in 27D, and 24 directions with 48 points each in 28D. These
are precisely the tail alphabets in the tetrad-head grammar.

For a fixed norm64 tail, norm128 points have squared-norm64 core vectors
with nonpositive pairwise products. A set of m unit vectors of this kind
in R24 has m≤48: for its Gram matrix G, the off-diagonal entries satisfy
Gij²≤-Gij. Positivity of the sum of all entries gives tr(G²)≤2m. Then
rank(G)≥tr(G)²/tr(G²)≥m/2, while rank(G)≤24.

Consequently both incumbent head layers saturate their entire tail alphabets.
No replacement using arbitrary real core bodies at these same tail vectors
can improve either total, even when all other incumbent points may be deleted.
This does not cover new tails, changed radii, or simultaneous new non-head layers.

## Global six-support midpoint grammar

`global_sixes.py` covers all C(24,6)=134,596 core supports, not just supports
already present in the source. A Golay/Steiner5 lookup prices the exact core
blockers. For supports contained in an octad, only that octad can conflict;
otherwise the six five-subsets identify all possible octads. All other core
shapes clear the checked largest-six-coordinate l1 bound. Each admitted
candidate's blocker list is independently rebuilt with direct int64 sums.

27D: 1,360,256 signed bodies pass the two-core-blocker gate, including support
families absent from the incumbent. Pairing with every signed tail root tests
16,323,072 candidates. The final two-blocker pool has 102 points. The numerical
MILP returns zero gain. This is a finite-pool outcome, not unrestricted optimality.
31D: 1,436,624 signed bodies paired with all 84 signed tail roots give
120,676,416 candidates. The final two-blocker pool has 750 new points;
numerical optimization returns zero gain on the 251,196-point source.
The factorized continuation reuses 92 completely scanned zero-output chunks
from the initial run; all admitted candidates are independently repriced
by direct integer coordinate sums. The scope remains the finite cutoff pool.

## Complete latitude-contact closure: the new result

`layer_contact_closure.py` enumerates negative-contact sums anchored on an
entire latitude, rather than the preceding campaign's 64-anchor samples.
In 31D the 6,592 old norm64-tail anchors yield 4,704,622 sum incidences and
1,243,881 distinct new points. The two-blocker pool has 52 candidates.
Selecting eight and deleting their four shared blockers gives 251,200.
The explicit orthogonal-block formula, complete integer coordinates, source
hashes, and independent audit are shipped in the portable package above.

A second complete pass from 251,200 returns no further gain. Complete norm64
latitude closures in 27D, 29D, 30D, and complete unpaired-anchor closures in
all 27–31D also return no gain in their two-blocker pools. These searches are
finite proposals, not global bounds. See their individual reports.

The complete tetrad grammars in 28D, 29D, 30D tested 4,080,384, 15,301,440,
and 42,844,032 candidates, respectively. The 28D case has the exact fiber
capacity proof above. The 29D cutoff-four MILP is unresolved. A strengthened
30D numerical model returns zero gain but supplies no exact solver proof.
The optional 31D tetrad sweep was stopped after 104 of 574 tail vectors;
its two saved partial pools do not constitute a complete tetrad sweep.

A final joint selection repriced the partial 31D tetrad pools and complete
six-support pool against the frozen 251,200-point witness. All 3,762 distinct
candidates have a blocker. An optimizer-free exact rational dual certificate
now proves that no exchange from this explicit pool can improve the size,
even with arbitrary old-point deletions. Its objective lower bound is
−428571624/1000000000 > −1, which excludes positive integral gain.
See [`joint31_exact/`](joint31_exact/) for the portable certificate and fresh
replay. This is a finite-pool closure, not a maximality or kissing-number upper bound.
