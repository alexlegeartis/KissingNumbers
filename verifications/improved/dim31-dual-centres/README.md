# Midpoints and dual E7 centres in dimension 31

Checkpoint: **26 October 2026**. Exact lower bound: **241,742**.
Previous method: [Felpix Studios, 241,079 points](https://github.com/Felpix-Studios/kissing-number/tree/b823640f4e98134b9d55b27288132ee954a8affa).
The retained local frame construction has 238,662 points.

196,560 − 43,264 + 86,528 + 1,792 + 126 = **241,742**.
338 midpoint blocks each gain 128 points. The additional 28 dual-centre
blocks each contribute 64 points, alongside the 126 pure E7 tails.
Midpoint rows have coordinates (x_0,...,x_29,sqrt(2)x_30).
Centre tails use this metric multiplied by sqrt(8/3). All norms are 32;
the product threshold is 16. `source/README31_241742.txt` supplies the full
construction, angle arguments, and attribution.

```bash
python final.py
python verify.py
```

The wrapper checks the complete SHA-256 manifest and coefficient ranges,
regenerates the Golay Leech shell, verifies owner membership and convexity
witnesses, and checks all midpoint, centre, and pure-tail pairs. Bounded
small integer dot products are represented exactly in binary64; radical
comparisons use signed integer squaring. No tolerance is used.
It independently materializes and compares every expanded coordinate array
with `source/points31_241742.npz`. Python and NumPy suffice.
Original receipts are retained; every replay uses a temporary directory.

Standalone scripts: `source/verify31.py` proves the angles and
`source/materialize31.py` expands the compact certificate.
