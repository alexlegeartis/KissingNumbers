# Explicit code layers and central sections

Checkpoint: **28 October 2026** (95–96D continuation). This package registers
eight constructions. The latest update improves its 95D and 96D certificates.

| Dimension | This construction | Mechanism |
|---|---:|---|
| 61 | **72,376,856** | binary latitude compression and sparse layers |
| 62 | **100,437,548** | two central sections of the 64D code construction |
| 63 | **182,537,988** | one central section |
| 66 | **338,302,340** | shortened weight-16 support code and lower layers |
| 67 | **351,420,806** | the same construction |
| 68 | **371,552,152** | the same construction |
| 95 | **7,092,331,276** | optimized central section of the enlarged 96D construction |
| 96 | **13,290,782,912** | transversal supports plus a mixed-pair support code |

These are constructive lower bounds. Cohn's reviewed October table stops at dimension
48, so it provides no comparison in this range. Historical comparisons remain
in [RESULTS.md](../../../RESULTS.md) and the supplied source reports. No
unrestricted optimum is asserted.

## Replay

The **28 October 49–72D continuation** adds a verified weight-eight bridge in
61D, parity-compatible sections of the augmented 64D parent in 62–63D, and
complete quadruple packings in 66–68D. Its gains over the preceding repository
certificates are 1,066,368, 66,080, 130,560, 10,240, 20,480 and 30,720 points.
[Full construction and proof scope](source/continuation49_72/README.md) gives the
coordinate formulas, inherited inputs, exact checks and discovery search.

```bash
python verify_continuation.py
```

This replay runs `source/continuation49_72/bridge61.py`, `signed_sections.py`
and `netto.py` in isolation. `bridge61.npz`, `section_shifts.json` and the
copied `parent64.json` are explicit proof inputs. Their JSON receipts are
recomputed, not treated as an optimization score. The optional
`search/seed.py` and `search/fill_correct.cpp` rebuild the bridge seed and
continue its search. The default `verify.py` includes all new checks.


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

In 96D, the current certificate supplies 24,645 additional weight-24 supports.
It starts from doubled QR48 supports and permits different pairings in
different four-cell rows, followed by neutral switches and additions. Each occupies at most twelve rows
of the existing 24-by-4 transversal layout, proving the cross-intersection
bound. The 95D construction is an exact central section, with sparse signs
ordered separately on each support. The source docstrings and JSON reports
give the count identities and geometric inequalities.

## Files

`final.py` derives the eight totals; `verify.py` replays their proofs.
`UPDATE_PROVENANCE.json` records research copies and documented changes to the
comparison fields. `data/SHA256SUMS.txt` covers the complete package except itself.

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
- `section95_search.py`, its JSON and NPZ reproduce the preceding 95D section.
- `support_pairing_certificate.py`, its JSON and NPZ independently rebuild
  the improved QR48 pairing and its best section. `support_pairing_best.txt`
  gives the coordinate permutation. The replay checks every QR48 word.

The older ERS, P48 cap and Gamma72 section packages remain valid retained
constructions. Their historical bounds are superseded by the values above.

## Optimizing the section directly

The 25 October CPU search weights each no-double support by its exact number
of signs surviving the selected two-row section. These multiplicities are
16,384, 8,192 and 6,144 according to whether the support occupies zero, one or
two selected rows. The objective includes the ordinary doubled supports too.
The preceding pairing certificate gives 233,595,904 support points in 95D
and 23,288 supports in its 96D parent. The continuation below raises both. The independent replay
checks the actual dense-code projection and all section populations; it does
not use an optimization score as a certificate.

## Mixed pairings in individual rows

The continuation replaces the uniform-partition restriction by an explicit
constant-weight support code. The first 28 October checkpoint had **24,631** supports, each of weight 24
and occupying at most twelve rows. Its verifier checked all **303,330,765**
distinct support pairs directly, with maximum intersection 12. It does not
trust a search graph or an optimization score. The improvement over the
preceding registered certificate (commit `a6ba338`) is **180,224** points in
95D and **294,912** in 96D. The previous bounds were **7,092,007,692** and
**13,290,258,624**, respectively. The support family replaces the earlier one;
its coordinates must not simply be appended to it.

For a short isolated replay, run `python verify95_96.py`. The ordinary
`verify.py` also includes this replay. `source/dim61_96/mixed_pair_certificate.py`,
its NPZ and JSON contain the complete support certificate and count proof.
The binary32 matrix products in this check are exact sums of at most 96
zeros and ones; no approximate angle comparison is used.

## Independent packed integer replay

`source/dim61_96/verify_support_bits.py` checks the same complete support
family using unsigned integer intersections and population counts, independently
of the coordinate matrix products. It checks all 303,330,765 unordered pairs,
weight 24, distinctness, at most twelve occupied rows and the exact section
histogram `[6999,0,11435,0,5463,0,704,0,44]`. No approximate arithmetic is used.
Its pinned receipt is `source/dim61_96/support_bits_verification.json`.
Both `verify95_96.py` and `verify.py` replay both representations in a temporary
copy and match the complete receipts (apart from elapsed time).

The current support family has 17,296 ordinary doubled QR48 supports and 7,349
compatible nonordinary supports. It contributes 403,783,680 points in 96D:

```text
8,589,934,592 + 4,294,967,296 + 403,783,680 + 2,097,152 + 192
    = 13,290,782,912.
```

The central-section terms are
`4,697,620,480 + 2,147,483,648 + 245,910,528 + 1,316,444 + 176
= 7,092,331,276`. These are constructive lower bounds, without an optimality
claim. `UPDATE_PROVENANCE.json` records the research source, exact certificate
hash and immediately preceding repository comparison.

The dense [96,33,24] code uses the nested concatenation described by
[Cheng's construction and its explicit code-table recipe](https://www.codetables.de/BKLC/BKLC.php?q=2&n=96&k=33).
The [binary code table](https://aeb.win.tue.nl/codes/binary-1.html) records the
16,384-word length-24 distance-six sign code; the replay enumerates that code
and the quaternary parent directly. Both source pages were reviewed on
27 October 2026. These code inputs are distinct from the new support family.
The manuscript and PDF retain their earlier bounds until the requested later
paper update.


## Fill/correct continuation on 28 October

The new family replaces 2,564 old supports by 2,582 supports, for a net gain
of 18. The final explicit certificate is the authority; no search bucket,
score or code-membership assumption is used by either independent verifier.

The finite search changes pairings inside four-cell rows and occasionally
moves a two-cell row to an empty row. A free candidate is inserted. Candidates
with exactly one blocker are grouped by that blocker; two mutually compatible
candidates can replace their shared blocker. Stored candidates are checked
again against the entire current family before a trade. Occasional one-for-one
moves explore configurations of the same size. Every accepted operation
preserves both the support count and the exact 95D section count, or increases
them. Integer population counts decide all conflicts.

This applies the fill/correct idea described by
[PackingStar](https://arxiv.org/html/2511.13391v4), whose Corrector removes
points before the next filling phase. The public repository supplies data
and verifiers, but no training implementation. This search does not reproduce
their reinforcement learning or their October 26–31D constructions.

The shipped starting checkpoint already contains two extra supports from
an earlier walk. From the package directory, reproduce its continuation:

```bash
python search/fill_correct.py --output /tmp/reproduced-supports.npz
```

The default seed is 282029, with 100,000,000 trials and neutral-period 256.
The continuation finds five free additions, eleven two-for-one trades and
89 one-for-one moves, ending at 24,631 supports. This produces the preceding checkpoint, which is the input to the further
continuation below; it is not the current certificate. NPZ compression bytes may
vary between environments. `search/receipt.json` records an ordered text
fingerprint and the search parameters. Searching requires C++17 and NumPy;
verification requires only NumPy. The search is not an exhaustion proof.

`python search/negative_controls.py` first accepts the current certificate,
then checks that both replays reject a duplicated support, a weight-23
support, and two supports with intersection 23. The recorded outcomes are
in `search/negative_controls.json`.

## Further support continuation, 28 October 2026

Continuing the verified 24,631-support family gives **24,645** supports.
Four free additions and ten two-for-one trades add fourteen supports, with
147 neutral exchanges. The bounds are **7,092,331,276** in 95D and
**13,290,782,912** in 96D, gains of **143,360** and **229,376** over commit
9d38a38. The 95D section count is optimized without lowering either bound.

Both complete verifiers checked all **303,675,690** distinct support pairs;
one uses physical incidence coordinates, the other packed integer
intersections. The sign-code distance enumeration and exact dense-code
section replay also passed. `verify95_96.py` runs both in an isolated copy.
`search/continuation_20261028.json` records the seed, input hash, search
parameters and result hash. The preceding `search/receipt.json` remains
historical; reproducing it supplies the starting checkpoint for this pass.
Search is finite and does not establish an optimum.
