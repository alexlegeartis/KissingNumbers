# A kissing configuration of 352320 points in ℝ³²

**τ₃₂ ≥ 352320**, obtained from a new binary constant-weight code showing **A(32, 8, 8) ≥ 1713**.

Qihong Lin (linqihong001@qq.com), 2026-10-22

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23527725.svg)](https://doi.org/10.5281/zenodo.23527725) — Zenodo archive, all versions: [doi:10.5281/zenodo.23527725](https://doi.org/10.5281/zenodo.23527725)

| | previous best | this repository |
|---|---|---|
| kissing number τ₃₂ (lower bound) | 347584 — Yang et al., [arXiv:2609.35051](https://arxiv.org/abs/2609.35051), as listed in [Cohn's table](https://cohn.mit.edu/kissing-numbers/) | **352320** (+4736) |
| A(32, 8, 8) (lower bound) | 1676 — same source, as listed in [Brouwer's table](https://aeb.win.tue.nl/codes/Andw.html) | **1713** (+37) |

Improved certificates, if any, will be added here as new releases (each archived on Zenodo).

New bounds in dimensions 33, 34, 35, 37, 38 and 39 (with the same method, plus an equatorial extension in dimension 38) are in a companion repository: [QH-Lin/kissing-numbers-33-39](https://github.com/QH-Lin/kissing-numbers-33-39).

## Construction

This is the three-layer construction of Edel, Rains and Sloane (*On kissing numbers in dimensions 32 to 128*, Electron. J. Combin. 5 (1998) R22). All vectors have squared norm 32.

1. **Dense layer** — 2¹⁷ = 131072 vectors ((−1)^{c_i})_{i} for every codeword c of the Cheng–Sloane [32, 17, 8] binary linear code (generator matrix in `gen32_17.txt`).
2. **Middle layer** — for each of the 1713 weight-8 supports S in `certificates/supports_1713.txt`, the 128 vectors equal to ±2 on S with an even number of minus signs and 0 elsewhere: 128 × 1713 = 219264 vectors.
3. **Root layer** — the 4·C(32,2) = 1984 vectors (±4, ±4, 0³⁰) (the D₃₂ roots, scaled).

Total: 131072 + 219264 + 1984 = **352320**. Any two distinct vectors have inner product ≤ 16, i.e. angle ≥ 60°.

The only new ingredient is the middle-layer support family: 1713 subsets of size 8 of {0, …, 31}, pairwise meeting in at most 4 points, i.e. a binary constant-weight code with n = 32, w = 8, d = 8.

### Why the inner products are at most 16

- dense–dense: 32 − 2·dist(c, c′) ≤ 16 because the code has minimum distance 8;
- dense–middle: at most 8 · 2 = 16; dense–root: at most 2 · 4 = 8;
- middle–middle, same support: two distinct even-sign patterns differ in at least 2 positions, so 4 · (8 − 2·2) = 16 at most;
- middle–middle, different supports: |S ∩ T| ≤ 4, so at most 4 · 4 = 16;
- middle–root: at most 2 · (2 · 4) = 16; root–root: at most 16.

So validity reduces to two facts, each checked separately below: the code has minimum distance ≥ 8, and the supports have weight 8 and pairwise intersections ≤ 4. The brute-force verifier does not rely on this argument; it checks every pair of the 352320 vectors directly.

## Verification

| program | what it checks | result for `supports_1713.txt` |
|---|---|---|
| `verify_code.py` | enumerates all 2¹⁷ codewords: dimension 17, minimum distance 8 | `[32,17,8]`, weight distribution 1, 908, 3328, 14784, 27392, 38246, 27392, 14784, 3328, 908, 1 (weights 0, 8, …, 32) |
| `verify_supports.py` | weight 8, all distinct, pairwise intersection ≤ 4 | 1713 supports, max intersection 4 |
| `verify_full.c` | builds all 352320 vectors and checks all 6.21 × 10¹⁰ pairs | max inner product 16, 0 violations (≈ 10 s on a 14-core Apple M-series laptop) |

To reproduce:

```bash
python3 verify_code.py gen32_17.txt
python3 verify_supports.py certificates/supports_1713.txt
cc -O3 -pthread -o verify_full verify_full.c        # on Apple silicon add -mcpu=native for NEON dot products
./verify_full gen32_17.txt certificates/supports_1713.txt
```

Expected last line of `verify_full`: `VALID kissing configuration of 352320 points in R^32`.

## Files

- `certificates/supports_1713.txt` — the 1713 supports (current best), one per line as a 32-character 0/1 string (leftmost character = coordinate 0). Older certificates in `certificates/` are kept for the record.
- `gen32_17.txt` — generator matrix of the Cheng–Sloane [32, 17, 8] code (Y. Cheng, N. J. A. Sloane, *Codes from symmetry groups, and a [32,17,8] code*, SIAM J. Discrete Math. 2 (1989) 28–37), transcribed from [codetables.de](https://www.codetables.de/BKLC/BKLC.php?q=2&n=32&k=17).
- `verify_code.py`, `verify_supports.py`, `verify_full.c` — the checks above.
- `SHA256SUMS` — hashes of all files.

SHA-256 of the current certificate:

```
539972293fbfc671477cea1805a5bf70a2208a4c1e60ec122347b5ef8a23c807  certificates/supports_1713.txt
```

## Record history

| date (UTC+8) | A(32,8,8) | τ₃₂ | certificate |
|---|---|---|---|
| 2026-10-21 14:47 | 1711 | 352064 | `certificates/supports_1711.txt` |
| 2026-10-21 14:58 | 1712 | 352192 | `certificates/supports_1712.txt` |
| 2026-10-21 15:29 | 1713 | 352320 | `certificates/supports_1713.txt` |

## Method (summary)

Finding the support family is a maximum independent set problem in the graph on all C(32, 8) = 10518300 eight-subsets of {0, …, 31}, two subsets adjacent when they share ≥ 5 points (degree 121264; the graph is never stored). Starting from the 1676-block family of Yang et al., we used an iterated local search in the style of Andrade–Resende–Werneck (insertions and (1,2)-swaps with a perturbation and acceptance rule), and then a ruin-and-recreate loop that deletes a large structured part of the current family and rebuilds it with the local search. The search code and a preprint describing it will be released later.

## Citation

Please cite the Zenodo archive: Qihong Lin, *A kissing configuration of 352320 points in R^32 and A(32,8,8) ≥ 1713: explicit certificates*, doi:[10.5281/zenodo.23527725](https://doi.org/10.5281/zenodo.23527725). See also `CITATION.cff`.

## License

[CC BY 4.0](LICENSE).
