# 37D: 536,440 by sharing signed-root deletions and adding a radical layer

Start with the [533,864-point QR witness](../dim37-38-qr-merge/).
Its dense and octad layers are unchanged. Replace 576 signed roots
with 720 stars, gaining 144, then add 2,432 points on 38 new signed supports. The total gain is **2,576**.

A star has one centre coordinate ±4 and four leaves ±2, so its
squared norm is 32. In each 16-point block the leaves run over all
eight even sign words, with both centre signs. First place nine
blocks on the special QR coordinate, one on each part of a partition
of the 36 ordinary coordinates into groups of four. This stage
replaces 144 signed roots with 144 stars. The next 36 blocks each
centre on an ordinary site, using the special coordinate as one
leaf. Their four signed roots on the special–centre pair are already
deleted. Each block adds 16 stars and uses at most 12 additional
ordinary signed roots. The actual recipe attains the net gain 144.

The [compact recipe](star_recipe.json) and verifiers check every
star–star pair, all affected signed root and octad geometry, and
exact dense layer norms and cross bounds. The second algorithm
builds the stars using even sign words instead of the Hadamard
rows, and verifies the signed-root cuts by incidence. The previous
certificate proves the retained interactions.

Run with Python 3 and NumPy:

```sh
python verify.py
python audit.py
```


## The seven-site layer

On each of the 38 six-coordinate ordinary supports in
`six_recipe.json`, place all 32 even-parity signings with
magnitudes 2. Independently set the special coordinate to
either +2 sqrt(2) or −2 sqrt(2). Each point has squared norm
32, and 64 points per support give 2,432 additional points.

The ordinary supports meet pairwise in at most two sites
and intersect any existing octad in at most four. Two
points with opposite special signs have product at most
−8+24=16; with equal signs, the overlap or even-sign parity
bounds their ordinary product by eight. No retained root
uses the special coordinate, so the seven-site layer incurs
no further root deletions. Star cross products are bounded
exactly as a+b sqrt(2). The 37D dense cross bound at special
magnitude m is 24+2m sqrt(2), at squared norm 144+m^2.

The main replay runs both stages. Independently,
`python audit_six.py` checks all 2,956,096 new same-layer
pairs and the exact cross-layer surd inequalities. Use
`python audit.py` for the first-stage sparse-root ledger.
