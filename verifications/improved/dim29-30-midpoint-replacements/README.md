# Midpoint replacements in dimensions 29 and 30

Checkpoint: **25 October 2026**.

| Dimension | Previous benchmark | This construction | Gain |
|---|---:|---:|---:|
| 29 | 210,495 | **216,568** | 6,073 |
| 30 | 222,073 | **224,280** | 2,207 |

The previous method is [K. Li's lifted Leech blocks and auxiliary layers](https://github.com/Felpix-Studios/kissing-number/blob/b823640f4e98134b9d55b27288132ee954a8affa/preprints/leech-lift-d28-d31/paper.pdf).
The earlier frame layers remain a retained construction in
[`../dim29-30-frame-layer/`](../dim29-30-frame-layer/).

At each chosen Leech residual, four deleted endpoints are replaced by eight
midpoint points with orthogonal tails of squared norm 8. All points have squared
norm 32 and products at most 16. Midpoint convexity proves compatibility with
every retained Leech point. Each uniform piece treats 32 residuals and gains 128
points. The 29D certificate uses 156 pieces and 40 pure D5 roots; 30D uses 216
pieces and 72 pure E6 roots. The inserted tails use D5 and D6 roots respectively.
The pure tail set can differ: Cauchy-Schwarz proves all mixed products.

```bash
python final.py
python verify.py
```

Python 3.9+ and NumPy suffice. The wrapper checks the complete SHA-256 manifest
and coordinate ranges, then runs independent proofs in a temporary copy.
The verifier regenerates the Leech shell, checks distinct deleted endpoints,
midpoint identities, exact norms, and every inserted/inserted pair. Float64 BLAS
performs only small integer operations with absolute partial sums at most 32;
there is no tolerance or rounding assumption. The source verifier also accepts
`--integer` for integer matrix multiplication.

`source/README.md` contains the proof and search provenance. The two compact
`certificate*.npz` files specify all replacements; `materialize.py` reconstructs
the full configurations as `a + b*sqrt(3)`. The expanded `points*.npz` files and
original verification receipts are included. Search optimality is not claimed.
`final.py` derives the counts; `data/SHA256SUMS.txt` covers all package files.
