# Addendum: larger shared-head supports in 22D

Use the same 49,656 frozen exterior rays obtained by releasing 240 rays incident to the head plane (ordinary coordinate t plus actual scalar height h). The definitions, benchmark and chamber walls in the companion README are unchanged. The exact theorem extends from at most five ordinary coordinates to **at most six arbitrary ordinary coordinates**. The only unit rays in that wider domain are still the original 240. No new bound.

## One parity relation still suffices

For an ordinary support J of size at most six, complete it to exactly six, and set S_0=J∪{t}. The projected sign code orthogonal is the shortening of D on S_0. Its minimum weight is at least five. If it contained two distinct nonzero words, the sum of the weights of those two and their sum would be at least fifteen but equal twice their union size, at most fourteen. Thus there is at most one parity relation, carried by a set H of size five, six or seven, and all sign assignments respecting it occur. D's weight-five words are the PG lines.

More generally the same at-most-one result holds for |S_0|≤8: the only possible weight triple with all weights ≥5 and even sum≤16 is 5,5,6. The two weight-five supports would be distinct PG lines and hence have union nine, contradiction. This observation is independent of the chamber replay.

## Uniform enlargement of the walls

Use ordinary magnitudes a_1≥a_2≥…≥0, z=|x_t|, Y=√3 h, and total T=z+Σa_i. The root, octad and one-coordinate outer walls remain
  a_1+a_2≤4; z+a_1+a_2+a_3≤8; a_1+Y≤8.
They follow exactly as in the companion proof: an octad containing t and the top three has a free sign coordinate (for |S_0|≤7), and among the four outer lines through the largest ordinary position avoiding t at most one is fully contained in the donor support. For |S_0|=8, the same octad wall follows because even if an octad is completely contained in S_0, there are **at least two** shortened octads through the given four positions and only one can equal S_0. Thus the familiar wall is valid through seven ordinary positions too.

If there is no parity relation or it has length other than five, no retained projective line is entirely contained in J. The largest ordinary position has some position in J on a different removed t-line (one t-line contains only four ordinary positions), so outer walls give a_1+a_|J|+Y≤8.
If H is a line containing t, take the largest ordinary coordinate inside it and the largest ordinary coordinate outside it. Their common line is retained and meets H only once. Since |J\H|≤3, that retained line is not contained in J, so use the sum of those two magnitudes plus Y≤8. If H is a line not containing t, retain its own two full sign banks and apply the two outer parity walls from the companion proof (sum_H+Y≤8 for p=0; sum_H−Y≤8 and sum_H−2min_H+Y≤8 for p=1).

For the middle walls write U=max_H a, v=min_H a, and B=max_{S_0\H} a when the outside set is nonempty; these include z wherever appropriate. For odd |H| the inside/outside parity walls are exactly those of the smaller line case, with the outside mark placed at B (ignore outside walls if no outside coordinate exists). For even |H| the complement does not switch parity, so the positive wing is dominant. If p=0 use
  T+2U−2v+Y≤16 and, if there is an outside mark, T+2B+Y≤16.
If p=1 use
  T+2U+Y≤16 and, if there is an outside mark, T+2B−2v+Y≤16.
For no relation use T+2max(a_1,z)+Y≤16. Penalties and zero coordinates are treated exactly as before.

Sorting z into the seven slots and choosing the outside positions of H (zero, one or two of seven), or no relation, gives 413 rational polytopes. The standalone fraction-free BigInt replay checks 2,495,493 row bases, including 216,317 feasible bases. Equality at norm 32 still occurs only at a released signed root or released outer line cap with its required parity; their magnitude patterns are the same as before, padded with an additional zero. Strict convexity therefore proves the stated theorem through six ordinary coordinates. No native computation or numerical solver was used.
