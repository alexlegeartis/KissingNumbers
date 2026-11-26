# 21D: every improvement of Li’s Golay line-label lift requires six-dimensional hole spans on both parities

The 21D record used for comparison is **31,560, Kevin Li (October 2026), not ours**. No larger configuration is claimed. This result concerns Li’s specified family of constructions only; it does not narrow the upper bound for kissing numbers.

Sources: Kevin Li, “Improved Kissing Numbers in Nineteen and Twenty-One Dimensions”, section 3, https://hexagonmath.org/2610.00248v1 ; independently reconstructed base and exhaustive pair check in https://github.com/alexlegeartis/KissingNumbers/tree/86cfd8463ac0adf5de10534bc26995c70f753f8f/research/compute19_21_20261104 . Existing independent-hole reductions, nonlinear examples and retained half-kernel obstructions: research/attack21_structured_20261104. Existing finite-family interval: 960 ≤ M ≤ 971 per component, research/upper21_certificate_20261104. Those prior results do not require holes to have affine dimension at least six on each separate parity. No claim of novelty beyond the checked repository/literature context.

## Geometry and finite lift

Let W be the 4096 signs obtained by deleting the last three coordinates of the extended binary Golay code. Its 21 weight-five words L are the lines of PG(2,4). Fix the 27,720 ray base: 840 signed two-entry roots with entries ±2 and 26,880 odd-signed unit-coordinate octads supported on the 210 surviving shortened weight-eight words. At a sign w∈W and label P∈L use the integer ray of magnitude 2 on P, 3 off P, and sign (-1)^w. Its squared norm is 164 and the normalized angle condition between these rays is precisely dot ≤ 82. Every candidate fits the fixed base: for the octads the maximum dot magnitude in the direction of matching signs is sum of magnitudes on the octad minus twice their minimum (Golay parity); line–octad intersections are 0, 2 or 4, giving maximum ≤18 < sqrt(328). For signed roots the maximum is at most 12. Distinct labels on one sign cannot coexist. The prior exact check on all label/sign differences proves: for different words at difference R∈L their two labels and R must be three distinct lines; all other sign differences cause no conflicts. The explicit dot check for all 4096 × 21² possibilities is in the linked reconstruction.

Thus the word graph has four interchangeable connected components (translations of span(L), size 1024); cross-component choices never conflict. The graph is bipartite by weight parity, with 512 words on each side, each an affine copy of V = F₂⁹. At most one ray per word. An improvement within this family requires a component with at least 961 labels, i.e. at most 63 unselected words (“holes”). Conversely a 961-word component copied to the four translates gives 27,720 + 4×961 = 31,564. A single improving component and three copies of Li’s 960 also suffice to obtain +1.

If x is selected with label P then x+P is a hole: otherwise R=P would be a conflict. If two adjacent selected words have the same label, they conflict. These two conditions are sufficient as well as necessary.

## Proposition (parity-flat barrier)

If the holes on either parity are contained in any affine subspace of dimension at most five, there are at least 64 holes in that component. No symmetry, affine assumption, or independence is imposed on the other side or on the assigned labels. The bound 64 is attained by Li’s affine-kernel construction. Therefore a component improving the published number must have affine hole span of dimension at least six on **each parity**.

Proof. By exchanging the sides and enlarging the flat, take an even affine 5-flat K containing all even holes E. Translate by an even vector to identify K with a linear 5-space. Choose any fixed line P₀ as odd origin. The quotient by K is F₂⁴ on each side; write a(P) for the quotient image of P+P₀ for P∈L. In these coordinates traversing the edge P adds a(P) to the syndrome. Write n_j = # {P : a(P)=j}; call j singleton if n_j=1. Set d=32−|E| and let F be the odd holes.

If two or more syndrome classes are missed by the directions, 64 odd words have no neighbor in K or E and must be holes. If exactly one is missed, 32 such holes are forced. There are at least 9 singleton syndromes (21 lines on 15 occupied classes). Each of the d points in K outside E used to be the unique K-neighbor of one odd point per singleton direction. These odd points are distinct, since an odd point with a singleton syndrome has just one K-neighbor. They are not in the missed class and are forced holes as well. Thus |E|+|F| ≥ 32−d+32+9d ≥64.

Suppose all 16 syndromes occur. After removing one copy of each from the 21 images, call the five remaining copies T. Every nonconstant affine halfspace in the quotient contains eight syndrome points. Its preimage captures 8+|T in halfspace| directions, hence a number between 8 and 13. The only possible numbers in this interval are 8, 9, 12, 13: the tiny Golay replay below checks all 512 linear characters of V, with complements giving both affine choices. Thus the length-five evaluation code of affine functionals on the five copies in T has only weights 0,1,4,5 and contains the all-ones word. Two distinct weight-one words would sum to a forbidden weight-two word. Every nonconstant evaluation is therefore the same singleton or its complement. Since affine functionals separate distinct quotient points, either the five copies coincide or four coincide and one lies elsewhere. The full profile of nontrivial multiplicities is consequently (6) or (5,2). In particular at most seven directions are in nonsingleton classes, and at least fourteen in singleton classes.

Consider x on the even side outside K; all 480 such x are selected. Suppose x has exactly one neighbor in F, reached along P. Then its own label must be P. Write s ≠ 0 for the syndrome of x and r=a(P). If r were singleton, choose a direction R with a(R)=s+r, using full syndrome coverage. R≠P. Then x+R is outside F, so selected, with syndrome r. Its only possible direction to K is P, hence its label must be P as well (if its unique neighbor is not in E it cannot be selected at all). This gives adjacent selected words both labeled P, a contradiction. Thus the sole-neighbor direction of x must belong to a nonsingleton class.

Let u count the outer even vertices with only one F-neighbor, and I all incidences from F to the 480 outer vertices. Each vertex has ≥1 neighbor in F, so 960 ≤ I+u. Full syndrome coverage gives each odd vertex a neighbor in K, whence I≤20|F|. Only seven directions can service a sole neighbor, so u≤7|F|. Consequently |F|≥ceil(960/27)=36. For d≤2 we obtain at least 66 holes in total. For d≥3, deleting a K-hole forces distinct additional odd holes in each of the at least 14 singleton-direction translates, so |F|≥14d and total holes ≥32+13d≥71. This finishes all cases. QED.

### Exact replay for the affine-hyperplane input

Run `python -B verify.py`. The complete 512-character histogram is checked from the original Golay generator.
```python
from collections import Counter
def span(generators):
    h = {0}
    for g in generators:
        h |= {w ^ g for w in h}
    return h
p = sum(1 << j for j in (0,1,5,6,7,9,11))
g = [(p << j) | (1 << 23) for j in range(12)]
words = {w & ((1 << 21)-1) for w in span(g)}
lines = sorted(w for w in words if w.bit_count() == 5)
assert len(words) == 4096 and len(lines) == 21 and len(span(lines)) == 1024
coord = {0: 0}
for l in lines:
    z = l ^ lines[0]
    if z not in coord:
        b = len(coord)
        coord.update({w ^ z: i ^ b for w, i in list(coord.items())})
assert len(coord) == 512
translated = [coord[l ^ lines[0]] for l in lines]
weights = Counter(sum((a & t).bit_count() % 2 for t in translated) for a in range(512))
assert weights == {0:1, 5:16, 8:130, 9:160, 12:120, 13:80, 16:5}
both = set(weights) | {21-w for w in weights}
assert set(range(8,14)) & both == {8,9,12,13}
```

## Constructive alternative for nonlinear hole sets

For fixed H=E∪F, give each h∈H a palette S_h of line directions. For adjacent holes require S_h and S_g disjoint. Cover each selected word w outside H by at least one palette: w+P∈H and P∈S_(w+P). This formulation is EXACT: assigning w one such label targets a hole and two adjacent equal labels P would translate to adjacent holes sharing P. Conversely a valid labeling defines these palettes by collecting labels into their target holes. A failed candidate is not a proof of a general bound.

Unlike the unrestricted 21,504 (word,label) graph, a specified 63-hole candidate uses at most 63×21 binary palette slots. Isolated holes may offer every direction without any side effects; independent H recovers the old sufficient domination construction. When the graph induced by H is a matching plus isolated vertices, a stronger exact shortcut applies: create one resource for each pair of incompatible slots and one free resource for each unpaired slot, connect each word to the resources of its possible slots. The slots’ conflict graph is a matching, since each slot’s target hole has degree ≤1. A matching covering every word is then equivalent to a valid selection. A failed Hall test applies only to that H; isolated user demands with no possible slots are failures immediately.

Asymmetric constructive sieve: fix E, compute the list D_y of directions from odd y to E. Odd y with empty list are forced holes. For each even x outside E and direction P the odd target x+P must be a hole if x takes P. Every neighbor y=x+R, R≠P, with D_y={P} must also be a hole. This is an inexpensive exact lower bound on the cost of a direction before trying an even assignment. If an even assignment is fixed (one label per selected even), an odd word is eligible iff no even targets it and at least one P in D_y does not appear at any selected even neighbor. All eligible odd words can then be taken simultaneously: odd–odd neighbors never conflict. This gives the exact best completion to any one-sided assignment.
