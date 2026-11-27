# 29D exact union comparison: no gain from the two existing witnesses

Lower86, 27 November 2026. Let A be [Kravatskiy's 232,157-point witness](https://github.com/alexlegeartis/KissingNumbers/tree/4e148069246da34e4326168e7e9f9dae745705ca/verifications/improved/dim29-even-grid-exchange), and B [Takhanov's 232,230-point v0.2 witness](https://github.com/k-nic/KissingNumbers/tree/77f0860bf93a745140e6e25481d3f60773b11441/version2). Their coordinate intersection is 172,246; the exclusive parts have 59,911 and 59,984 points. All coordinates use exact squared norm 128 and dot cap 64.

The exposed bipartite conflict graph has an explicit matching of size **59,911**, saturating the smaller side. Every matching edge has raw dot at least 80, and the certificate uses each second-side vertex only once. Therefore any kissing subset of the union can contain at most one endpoint per matched edge. There are 172,246 common points, 73 unmatched vertices and at most one possible choice per edge, giving the upper bound 172246+73+59911=232230. The witness B attains it. No claim is made about another geometry or a perturbation of either input.

The portable replay `verify-union.mjs` takes these five inputs in order:
1. [Takhanov v0.2 d29_certificate.npz](https://github.com/k-nic/KissingNumbers/blob/77f0860bf93a745140e6e25481d3f60773b11441/version2/certificates/data/d29_certificate.npz), SHA256 b27ec26a1d186ed6e9b6885921abf3b2e6bbde132206fde2689dfdaa7e39f824.
2. [Takhanov v0.1 d29_alt_certificate.npz](https://github.com/k-nic/KissingNumbers/blob/dadcb2940f7c3c3bd478b0543931d3eef9f4bd65/version1/certificates/data/d29_alt_certificate.npz), SHA256 ded2da5a8141af3606d375df464f2ce3bfa0c45bb02204610d94d19ccdcf8cd5.
3. [Own +1 exchange from the v0.1 alternant](https://github.com/alexlegeartis/KissingNumbers/blob/4e148069246da34e4326168e7e9f9dae745705ca/verifications/improved/dim29-takhanov-auxiliary-exchange/data/exchange29.npz), SHA256 130c0c144b014e1afa2284242718918fa7ccc2a9af8dd026a10d1a6a5aa9bd3e.
4. [Own +3 even-grid trade](https://github.com/alexlegeartis/KissingNumbers/blob/4e148069246da34e4326168e7e9f9dae745705ca/verifications/improved/dim29-even-grid-exchange/trade.npz), SHA256 d85b5cdf65eb128ad55104711cfc57621554de1d31005ad181757b46b2f3f9cd.
5. `matching-edits.json`, a literal replay stream of only 2,675 edits to the 57,385 same-body greedy edges, not a maximum matching solver.

It also accepts the separately published own 232157-coordinate NPZ as an optional sixth input and compares the literal published row hash and the reconstructed canonical coordinate set. The intermediate own 232154 canonical hash is f62df9eb2b3cdf7442828722aeb8b7e66441ffc6c1a16043429f88bec0904a3a (matches its published prior proof), and the final own 232157 canonical hash is 85d48a37bccd4aa5950c3b5279b9a673511a71ebc1eef7f0df9647fa1d32495d. The official B canonical hash matches Takhanov's manifest: 3675c465c219c0e542f881356774f85dc9d67cb7ad0f83cfb805d0e808138233. The worker checked all 59,911 exact integer edges and source hashes hosted. An independent comparison to the original full prior coordinate file reproduced both its row-order and canonical set hashes. No new record.

## Distinct radial scope on four exact old bodies

Take the body projections of the last four old-only vertices left after all short alternating paths (their own old tail profiles were three axial rays and one D5 midpoint). At each body b allow a nonnegative common body scale u and a nonnegative tail coefficient v along any signed equal-amplitude direction t of support size 1, 2, 4 or 5. There are 10+40+80+32=162 tail lines per b, including both orientations as rays.

For each of the 172 distinct actual tail vectors q in Takhanov's witness, compute P(q)=max b·p over all points (p,q) in the full exact reconstructed B, with Leech tail zero inclusively. Then the whole exterior imposes the rational polygon u P(q)+v (t·q)≤64. Intersect with 0≤u,v≤16: this does not lose a candidate of squared norm 128, since both component norms squared are positive integers. Perform exact homogeneous polygon clipping, retaining each old vertex in the halfplane and each literal rational edge crossing. Since the norm-square form is positive definite, its maximum is at a vertex. The maxima over the 162 choices are 80,80,80,512/5 for the four bodies, strictly below 128 in all 648 exact tests. Thus common radius variation cannot salvage these particular old-only donors at these directions. This does not rule out any arbitrary tail/contact move on other bodies or a profitable replacement outside the input union.

The second replay is `verify-radial.mjs`, with only the same pinned v0.2 compact NPZ as input.

## Scoped radial diagnostic

`verify-radial.mjs` also checks four pre-augmentation donor bodies against the 232,230-point witness. For each body it tests all 162 signed equal-amplitude directions of support sizes 1, 2, 4 or 5, leaving the body and tail amplitudes free and nonnegative. The exact optimum squared length is at most 80, 80, 80 and 512/5 for the four bodies, respectively, below the required 128. This does not exclude other bodies or arbitrary tail directions.
