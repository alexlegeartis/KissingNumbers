# Joint latitude families in dimension 38

Checkpoint: **26 October 2026**. Exact lower bound: **610,394**.
Previous method: [Qihong Lin's Leech cap construction](https://github.com/QH-Lin/kissing-numbers-33-39/tree/0c4e2bd486061b3fd55e128ac3497bfb12b73f07), **605,322** points.
The retained local rational latitude construction has 605,346 points.

589,680 caps + 13,710 equatorial points + 1,568 retained flat axis points
+ 28 coordinate axis points + 4 × 1,352 new latitude families = **610,394**.
Body lines may repeat. Every pair of new families satisfies the joint
integer inequality |U·U'| + 4|Z·Z'| ≤ 64. The signed diagonal orthogonal
axis transformation and all old/new angle types are checked exactly.

```bash
python final.py
python verify.py
```

The complete SHA-256 manifest covers the independent integer verifier,
certificate, original receipt, and all inherited input data and attribution.
The proof in `source/proof.txt` establishes norms, pair inequalities and
distinctness. Python and NumPy suffice; no optimizer, cache, download or
tolerance is needed. The verifier rebuilds the Leech shell and replays in
a temporary directory.
