# The 32-donor 26D variation is capped at 19 chamber rows

This is for a specific exact exchange against [Takhanov v0.2](https://github.com/k-nic/KissingNumbers/tree/77f0860bf93a745140e6e25481d3f60773b11441/version2) at 201,212. It does not decide the alternative [24-donor gate](four_octahedra_21.md). Coordinates in each chamber are hot u, three regular coordinates r, old tail v, and an alternate tail. The chart orders, including the two tails, are (0,3,1,2,9,8), (0,3,5,7,8,9), (4,6,1,2,8,9), (4,6,5,7,9,8). In each place take eight donor rows (±7, ±5, ±6, ±6, 7, 1)/14, with hot sign independent and even regular signs. Integer replay uses the [pinned upstream chart file](https://github.com/k-nic/KissingNumbers/blob/77f0860bf93a745140e6e25481d3f60773b11441/version2/certificates/data/d26_blocks.json).

These 32 donors are mutually legal; numerator dot over all pairs is at most 98. Axis and exterior caps follow from componentwise 7 and regular head sum plus twice the alternate tail equal to 26≤28. The only old skeleton payments are the four signed rows of tetrad [0,4,8,9] with both tails positive. With another original five-coordinate support the sum of donor numerator magnitudes is at most 13<14, so the caps suffice for cross-chamber interactions. Distinct clipped five-coordinate chambers and the old tetrads intersect any such chamber in at most two coordinates. The frozen count is 200864+20+224−4+32=201136.

Each rebuilt chamber is a clipped unit code in five coordinates (u,r₁,r₂,r₃,v), with all magnitude caps 1/2 and extra cutoff
`7|u|+7v+Ψ(r)≤7`,
where Ψ is the maximum over even signings of 5r₁+6r₂+6r₃. At most one code site fits each sign orthant (aligned dots≥3/4), and any occupied orthant adjacent to another occupied one across coordinate i has |xᵢ|≥1/(2√2) by the [directed rule](takhanov26_2513_23.md).

## What a 20-row code would have to look like

In the northern half v≥0 all eight orthants of even regular parity are forbidden: the full score is ≥5 times the sum of magnitudes, which is at least 2. In an odd regular orthant the maximum is weighted sum minus twice its minimum, and is at least \(L\sqrt{1-u^2-v^2}\), where \(L=\sqrt{450/43}\) follows from the maximum of the three weighted magnitudes. If both |u| and v are ≥1/4, the resulting score lower bound \(7(|u|+v)+L\sqrt{1-u^2-v^2}\) is increasing in either variable on their clipped square. Norm requires u²+v²≥1/4, so its minimum is over the lower arc, at least \(7(1+\sqrt3)/4+L\sqrt3/2>7\). Thus a northern odd row has |u|<1/4 or v<1/4, prohibiting its respective adjacent orthant. If both northern hot signs are occupied for the same regular signs, neither |u| can be small, so both corresponding southern mates must be absent. Consequently each odd regular sign word has at most three of its four hot/tail sign states. This also shows a 20-row code would occupy all sixteen southern sign orthants, plus one northern site for each of the four odd regular sign patterns, with |u|<1/4 and v≥1/4.

## The northern sites would force all sixteen southern tails too high

For any required northern site put \(W=u^2+v^2-1/4\). Then 0≤W≤1/16 and all three regular magnitudes are ≥√3/4. The sum of their deficits from 1/2 is at most \(W/(1/2+\sqrt3/4)\), since their squared magnitudes sum to 3/4−W. The smallest weighted regular magnitude is at most 5/2, so
\[
\Psi(r)\ge 7/2 - 6W/(1/2+\sqrt3/4).
\]
Meanwhile \(|u|+v-1/2\ge\sqrt{1/4+W}-1/2\ge W/(1/2+\sqrt5/4)\). Since \(7/(1/2+\sqrt5/4)>6/(1/2+\sqrt3/4)\), the cutoff forces W=0, |u|+v=1/2 and hence u=0,v=1/2, all regular magnitudes 1/2.

With u=0 such a northern site is geometrically sign-aligned on u with **either** of its two southern mates with the same regular sign word. Applying the directed rule on their sole sign difference v gives t=−v≥1/(2√2) for both mates. These are all eight odd southern sites. The other eight are even; for them every non-tail magnitude is ≥√(1/4−t²), and the full score is ≥24√(1/4−t²)−7t, above 7 when t≤3/10. Hence all sixteen southern sites have 3/10<t≤1/2.

The sixteen occupy an entire sign 4-cube on (u,r₁,r₂,r₃), so their magnitudes in all those coordinates are ≥1/(2√2). Apply the [4-cube edge-energy argument](takhanov26_48_donors.md) with all vertices high: with deficits d_i=1/2−|x_i|, one has t²=Σ(d_i−d_i²). The sign-flipping edge side summed over the cube is ≥2Σt². The four-coordinate variation is <0.516Σt², and tail variation is ≤(8/17)Σt² because \((t-t')^2\le(2/17)(t^2+t'^2)\) on [3/10,1/2]. This is a contradiction. Thus a rebuilt chamber has at most 19, and even four independently rebuilt chambers give at most 201136+4·19=201212. No new numerical bound is claimed.

Finite replay: `python3 verify_32_donors.py path/to/d26_blocks.json`.
