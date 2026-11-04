# Certificates for new kissing-number lower bounds in dimensions 25, 26, 27, 29 and 31

This repository contains the exact certificates and verification scripts for

> R. Takhanov and S. Yun, *New lower bounds for kissing numbers in dimensions 25, 26, 27, 29 and 31*, 2026. arXiv:XXXX.XXXXX <!-- TODO: arXiv ID -->

| d  | previous record | new bound  | gain | construction | folder |
|----|-----------------|------------|------|--------------|--------|
| 25 | 197580 [YZW+]   | **197583** | +3   | new heads + coordinated motion of 552 bulk points | `d25/` |
| 26 | 199806 [K]      | **199868** | +62  | deformed second cap layer + 22 free heads | `d26/` |
| 27 | 201567 [YZW+]   | **201589** | +22  | larger second cap layer (graph colouring) | `d27/` |
| 29 | 209968 [K]      | **209970** | +2   | rotated lifted directions, 34-point auxiliary block | `d29/` |
| 31 | 238662 [K]      | **238666** | +4   | rotated lifted directions and rotated E7 auxiliary block | `d31/` |

Every bound is certified in exact arithmetic. No floating-point tolerance is used in any acceptance decision.

## Quick start

Requirements: Python ≥ 3.9 and NumPy. Nothing else is needed and nothing is downloaded. The scripts were tested with Python 3.12 and NumPy 2.4.

Each verifier can be run from any working directory:

```bash
python d25/d25/verify25.py        # ~5 s
python d26/verify26_free.py       # ~20 s
python d27/verify27_layer2.py     # ~20 s
python d29/verify29.py            # ~2 s
python d31/verify31.py            # ~2 s
```

Times are for one CPU core. Peak memory is under 1 GB.

A successful run ends with:

| script | success line |
|--------|--------------|
| `verify25.py` | `"status": "EXACT VERIFICATION PASSED"` and `"points": 197583` |
| `verify26_free.py` | `ALL CHECKS PASS   tau_26 >= 199868` |
| `verify27_layer2.py` | `ALL CHECKS PASS   tau_27 >= 201589` |
| `verify29.py` | `*** ALL CHECKS PASS ***   tau_29 >= 209970` |
| `verify31.py` | `*** ALL CHECKS PASS ***   tau_31 >= 238666` |

Every verifier exits with a non-zero status on failure. All five can be run at once:

```bash
for s in d25/d25/verify25.py d26/verify26_free.py d27/verify27_layer2.py d29/verify29.py d31/verify31.py; do
  python "$s" > /dev/null && echo "PASS  $s" || echo "FAIL  $s"
done
```

The output of the run made when the package was assembled is stored next to each script as `verify*.log`.

### Optional floating-point cross-checks

Each folder also has a `fullcheck*.py`. It builds all coordinates in float64 and computes every pairwise inner product. These scripts are an independent safety net, not part of the proof. They take from about a minute up to about 7 minutes (`d31`).

```bash
OPENBLAS_NUM_THREADS=1 python d25/d25/fullcheck25.py
python d25/d25/fullcheck25.py --export configuration_197583.npz   # explicit (197583, 25) coordinate array
```

## Repository layout

```
d25/
  README.md  SHA256SUMS.txt  requirements.txt
  d25/
    cert25_197582.json       precursor: 3 class shifts + owner lists, 48 local heads (197582 points)
    cert25_197583.json       extension: layer motion, added point, 3 replacement heads
    leech25.py               integer Leech minimal shell from the Golay code
    verify25.py              exact verifier
    fullcheck25.py           float cross-check / coordinate export
    verify25.log  verification.json  fullcheck25.log
d26/
  cert26_199868.json         22 free heads, 320 second-layer heads, 4 deleted owners
  data/heads26_Y.npy         first layer (1458 heads), from [K]
  data/heads26_side.npy
  verify26_free.py  fullcheck26.py  README.md  *.log
d27/
  heads27_layer2_u.npy       333 second-layer owners (new certificate)
  heads27_layer2_line.npy    their lines 0/1/2 = e1/e2/e3
  data/heads27_Y.npy         first layer (2258 heads, B. Lindow), via [K]
  data/heads27_side.npy
  verify27_layer2.py  fullcheck27.py  README.md  *.log
d29/
  axis34_rotated_caps.json   rational rotation R and 34 rational auxiliary vectors (new certificate)
  data/owners29.npy  data/bounds29.npy  data/geom29.json   fixed data from [K]
  verify29.py  fullcheck29.py  README.md  verify29.log
d31/
  axis122_d31.json           Cayley data for Q, R and the 122 kept roots (new certificate)
  data/owners31.npy  data/bounds31.npy  data/geom31.json   fixed data from [K]
  verify31.py  fullcheck31.py  README.md  *.log
```

Each folder's own `README.md` documents the file formats in detail.

### Checksums of the certificate files

```
238f68a77e2d5373c009c0eaad95468b9f7dc627a04671afe7b537d7b3b5193e  d25/d25/cert25_197582.json
bd3a684ef42342a573fcab7c7a6cbebfe87f99a76a408a0ee6f90412bb9b4a87  d25/d25/cert25_197583.json
2a11240379242abaf54e6e19107ef7598e8c7b0f301abe8057fa63b7dc16b9c8  d26/cert26_199868.json
91876a719499a9002b68c6f94b76fa819753262e854530d19a67d0af21a8d0b2  d27/heads27_layer2_u.npy
d16172cc8515649acf19ae5da83b9ec8cab55da89feab5abd1cfc777dc54a820  d27/heads27_layer2_line.npy
af9df2dd2d373cf6b4dd1d819e0139813ad52667e74fb64b9c885cf5b80bae8e  d29/axis34_rotated_caps.json
a5850d521b246fe1d6cb5cf7380d622943ac98ecd24cb54c178e06ef05a7c993  d31/axis122_d31.json
```

`d25/SHA256SUMS.txt` additionally covers every file in the `d25` package. `cert25_197583.json` records the SHA-256 of its precursor, and `verify25.py` checks it.

## Common conventions

**Leech lattice.** All five verifiers rebuild the 196560 minimal vectors in the same coordinate model. They start from the cyclic Golay code of length 23 with generator polynomial exponents `[0, 2, 4, 5, 6, 10, 11]`, extended by a parity coordinate. The lattice is scaled so that minimal vectors have squared norm 32, with shapes (±4², 0²²), (±2⁸, 0¹⁶) and (∓3, ±1²³). Owners and heads are always identified by their integer coordinate rows, never by shell indices.

**Kissing condition.** A configuration is a set of unit vectors with pairwise inner products at most 1/2. Some scripts work in units where every point has squared norm 32 and the condition reads ⟨P, Q⟩ ≤ 16.

**Exact arithmetic.** Inner products on the Leech side are integers. Quantities in ℚ(√2, √3) are stored as integer quadruples (a, b, c, d) representing (a + b√2 + c√3 + d√6)/m. Their signs are decided exactly or with rigorous integer bounds for √2, √3, √6 at precision 10⁻⁴⁰ to 10⁻⁸⁰. In `d26`, terms involving √N for an integer N are bounded the same way.

## The five constructions

The paper gives the full definitions and the lemmas that reduce the kissing condition to the inequalities checked here. "Bulk" or "equator" means the normalized Leech minimal vectors (z/√32, 0) that are not removed.

### d = 25: τ₂₅ ≥ 197583

Heads are vectors of squared norm 3/4 in ℝ²⁴. Each head removes one owner from the bulk and carries the two caps (x, ±1/2). There are two kinds of heads:

- **971 class heads** x = (U + tV)/√32, with t = (3 − √3)/6 and three shifts V₁, V₂, V₃ of squared norm 48. The owner sets have sizes 552, 286 and 133.
- **48 local heads** x = (A + √6 B)/(q√32), with AᵀB = 0 and |A|² + 6|B|² = 24q².

Together with the poles (0, ±1) and the point (−V₁/√48, 0), these give the 197582-point precursor.

The extension follows the coordinated motion of [YZW+, §3]. The 552 bulk points b with V₁ᵀz = 24 move to (b + δp, h), with h = (√6 − 2)/8 and δ = √6/4 − √3/4 − √2/8. These are the antipodes of the class-1 owners. The final three local heads are replaced (same owners), and the point (V₁/8, −1/2) is added.

| block | count |
|-------|-------|
| fixed bulk | 194989 |
| moved bulk | 552 |
| caps (2 × 1019 heads) | 2038 |
| poles | 2 |
| (−V₁/√48, 0) and (V₁/8, −1/2) | 2 |
| **total** | **197583** |

### d = 26: τ₂₆ ≥ 199868

This construction modifies the configuration of [K]. The first layer consists of points (Y/(3√32), t/√3). Here Y has squared norm 192, and t runs over one of two triangles: side 0 at 30°, 150°, 270°, or side 1 at 90°, 210°, 330°. The auxiliary hexagon has vertices at 0°, 60°, …, 300°.

The 22 new **free heads** are Y = 2v with v of squared norm 48. Since ⟨Y, z⟩ ≤ 48 on the whole shell, they remove no bulk point. To make room for them, the second layer is re-optimized: 4 of its heads are deleted, and their owners return to the bulk.

Each second-layer head is stored as an integer vector x = (p, q) ∈ ℤ²⁴ × ℤ². It contributes the points (p, ±q)/|x|; no norm is assumed.

| block | count |
|-------|-------|
| equator | 194782 |
| first layer (3 × (1458 + 22)) | 4440 |
| second layer (2 × 320) | 640 |
| hexagon | 6 |
| **total** | **199868** |

### d = 27: τ₂₇ ≥ 201589

The first layer of 2258 heads (2090 with one owner each, 168 free) and the 12-point auxiliary configuration are kept unchanged. The first layer is B. Lindow's, as distributed in [K]. The auxiliary configuration consists of the cuboctahedral roots rotated by 45° about e₃.

The second layer consists of points (√3u/(2√32), ±e_ℓ/2). It is admissible if and only if

- ⟨Y, u⟩ ≤ 32 for every first-layer head Y, and
- ⟨u, u′⟩ ≤ 8 for owners on the same line.

So choosing it is a proper 3-colouring problem on minimal vectors. The certificate takes 333 owners: 113, 112 and 108 on e₁, e₂, e₃.

| block | count |
|-------|-------|
| equator | 194137 |
| first layer (3 × 2258) | 6774 |
| second layer (2 × 333) | 666 |
| auxiliary | 12 |
| **total** | **201589** |

The two `heads27_layer2_*.npy` files use the format of the KissingNumbers package [K]. They can therefore also be checked with that package's `verify27.py`, which requires SymPy.

### d = 29: τ₂₉ ≥ 209970

The fixed data come from Kravatskiy's frame-layer configuration [K]:

- 14 disjoint classes of 248 type-B owner lines;
- the normalized D₅ roots, split into 12 zero-sum triangles and 2 antipodal pairs;
- the 480-point frame layer (f/√2, ±e_j/√2).

The certificate applies a rational rotation R ∈ SO(5) to the additional coordinates of all lifted vectors. Jointly with it, a new auxiliary block A of 34 rational unit vectors replaces the original 32. Taking R = I and the 32-point auxiliary block stored in `geom29.json` recovers Kravatskiy's 209968-point configuration.

| block | count |
|-------|-------|
| bulk | 189616 |
| lifted (√(2/3)u, √(1/3)Rt) | 19840 |
| auxiliary (0, a) | 34 |
| frame layer | 480 |
| **total** | **209970** |

### d = 31: τ₃₁ ≥ 238666

The fixed data again come from [K]:

- 42 disjoint classes holding 10328 type-B owner lines;
- the normalized E₇ roots, split into 42 zero-sum triangles;
- the frame W and the 672-point frame layer.

The lifted directions are rotated by R. The auxiliary block consists of the images under Q of 122 E₇ roots. Both rotations are exact rational Cayley transforms, (DI − K)(DI + K)⁻¹ with D = 10¹² and integer skew-symmetric K. Exactly 122 roots are admissible for the certified pair (Q, R), against 118 in the original configuration.

| block | count |
|-------|-------|
| bulk | 175904 |
| lifted | 61968 |
| auxiliary (0, Qe) | 122 |
| frame layer | 672 |
| **total** | **238666** |

## What the verifiers check, and what they assume

Each verifier rebuilds the Leech shell and the fixed data, then checks the certificate against the inequalities of the corresponding lemma in the paper. This covers norms, owner uniqueness and distinctness, and every head–head, head–bulk, layer and auxiliary condition. It then recomputes the point count.

The verifiers rely on two inputs:

1. **The rebuilt shell is the set of minimal vectors of the Leech lattice.** This is the standard Golay-code construction. In particular, distinct bulk vectors have inner product at most 16; `verify27_layer2.py` also checks the Leech inner-product distribution against a fixed vector.
2. **The reduction lemmas of the paper.** They show that all pairs not listed explicitly satisfy the kissing bound automatically. Examples are bulk–bulk pairs and pairs whose inner products are preserved by an orthogonal map.

The `fullcheck*.py` scripts do not use the lemmas. They test every pair numerically and confirm the same conclusion in floating point.

You can also check the fixed data against the upstream KissingNumbers repository instead of the copies shipped here:

```bash
git clone https://github.com/alexlegeartis/KissingNumbers
python d29/verify29.py      --repo ./KissingNumbers
python d31/verify31.py      --repo ./KissingNumbers
python d26/verify26_free.py --repo ./KissingNumbers
python d27/verify27_layer2.py --repo ./KissingNumbers
```

## How the constructions were found

All new constructions were generated with AI systems: Claude Code with Claude Opus 5.5, and GPT-6 Astra. For d = 29 and d = 31, the systems were asked to rotate the additional coordinates of the lifted vectors and fit more points into the auxiliary block. That idea was first used in d = 31 in [TY]. For d = 25, 26 and 27, they were asked to improve the existing records. Correctness rests entirely on the exact verifiers above, not on the search procedure.

## References

- **[K]** A. Kravatskiy, *New lower bounds for kissing numbers in dimensions 18 through 96*, preprint, September 2026; verification packages v1.1.1 at <https://github.com/alexlegeartis/KissingNumbers>.
- **[YZW+]** S. Yang, R. Zhao, J. Wu, et al., *Large-scale autonomous discovery of kissing number constructions*, arXiv:2609.35051, 2026.
- **[TY]** R. Takhanov and S. Yun, *New lower bounds for kissing numbers in dimensions 25–31*, arXiv:2609.21591, 2026.
- **[PS]** C. Ma et al., *Finding kissing numbers with game-theoretic reinforcement learning* (PackingStar), arXiv:2511.13391.
- **[CJKT]** H. Cohn, Y. Jiao, A. Kumar, S. Torquato, *Rigidity of spherical codes*, Geom. Topol. 15 (2011), 2235–2273.

## Citation

```bibtex
@misc{TakhanovYun2026kissing2531,
  author = {Rustem Takhanov and Stanislav Yun},
  title  = {New lower bounds for kissing numbers in dimensions 25, 26, 27, 29 and 31},
  year   = {2026},
  note   = {arXiv:XXXX.XXXXX}  % TODO
}
```

## License

TODO: choose a license (for example MIT for code and CC BY 4.0 for data). The fixed data in `d26/data`, `d27/data`, `d29/data` and `d31/data` are redistributed from [K], which is MIT-licensed; keep that notice.

## Contact

Rustem Takhanov, Mathematics Department, Nazarbayev University, Astana, Kazakhstan. <rustem.takhanov@nu.edu.kz>
