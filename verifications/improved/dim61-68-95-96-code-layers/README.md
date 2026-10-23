# Explicit code layers and central sections

Checkpoint: **24 October 2026**. These constructions replace the earlier
repository values in six dimensions and add dimensions 66 and 67,
registering eight constructions in total.

| Dimension | This construction | Mechanism |
|---|---:|---|
| 61 | **71,310,488** | binary latitude compression and sparse layers |
| 62 | **100,371,468** | two central sections of the 64D code construction |
| 63 | **182,407,428** | one central section |
| 66 | **338,292,100** | shortened weight-16 support code and lower layers |
| 67 | **351,400,326** | the same construction |
| 68 | **371,521,432** | the same construction |
| 95 | **7,077,231,372** | optimized central section of the enlarged 96D construction |
| 96 | **13,264,027,840** | transversal supports plus doubled QR48 supports |

These are constructive lower bounds. Cohn's current table stops at dimension
48, so it provides no comparison in this range. Historical comparisons remain
in [RESULTS.md](../../../RESULTS.md) and the supplied source reports. No
unrestricted optimum is asserted.

## Replay

```bash
python final.py
python verify.py --full-codes
```

Python requires NumPy; `--full-codes` additionally requires `g++`. The default
`verify.py` replays the binary compression, the [64,28,16] enumeration,
the sparse layer constructions and every central-section count. It reads the
recorded [73,36,16] distance certificate. `--full-codes` freshly checks all
**68,719,476,736** words of that code, comparing its complete weight histogram
with the supplied report. Both modes work in a temporary copy and preserve
the installed artifacts. `common/checkpoint.py` supplies the manifest check
and isolated replay helpers.

## Mathematical scope

The 61D map compresses the explicit [62,26,16] code into 61 coordinates with
two head magnitudes. Exact rational and quadratic-radical comparisons check
every permitted fibre/distance case. A quaternary [15,6,8] support code,
binary [15,10,4] sign code and all weight-two points add compatible layers.

For 62D and 63D, weight-eight dual words select central hyperplanes of the
64D construction. Exact projection counts and support-dependent sign orders
give the surviving populations. The 66–68D constructions use explicit
shortenings of the [73,36,16] weight-16 shell, extended Hamming signs and
affine weight-four support packings. Code distances bound all cross-support
angles; sparse/dense cross bounds are checked exactly.

In 96D, 23,012 additional weight-24 supports come from doubled QR48 supports
and two compatible doubling partitions. Each occupies at most twelve rows
of the existing 24-by-4 transversal layout, proving the cross-intersection
bound. The 95D construction is an exact central section, with sparse signs
ordered separately on each support. The source docstrings and JSON reports
give the count identities and geometric inequalities.

## Files

`final.py` derives the eight totals; `verify.py` replays their proofs.
`UPDATE_PROVENANCE.json` records byte-preserving copies from the research
workspace. `data/SHA256SUMS.txt` covers the complete package except itself.

Under `source/dim49plus/`:

- `binary_62_code.py` rebuilds the binary generators and enumerates weights.
- `latitude_compression.py`, `latitude_compression.json` give the 61D map,
  exact geometry, negative control and recorded counts.

Under `source/dim61_96/`:

- `weighted_layers.py`, its JSON and NPZ give the explicit binary/quaternary
  codes, nested concatenation proof and initial sparse layers.
- `source_sign23.json` preserves the sign-code input and provenance.
- `dual_sections.py`, its JSON and NPZ give the two central sections.
- `code64_supports.npz` contains the exhibited weight-16 supports; membership
  is checked against the rebuilt generator.
- `expanded_support_layers.py`, its JSON and `shortening_support_search.json`
  give the selected 66–68D shortenings and exact support counts.
- `cyclic73.generator.txt`, `cyclic73.weight16.bin`, `cyclic73.json` contain
  the binary generator, weight-16 supports and complete weight histogram.
- `enumerate_binary.cpp` independently enumerates that binary code.
- `support_doubling.py`, its JSON and NPZ construct the additional 96D supports.
- `section95_search.py`, its JSON and NPZ reproduce the strongest 95D section.

The older ERS, P48 cap and Gamma72 section packages remain valid retained
constructions. Their historical bounds are superseded by the values above.
