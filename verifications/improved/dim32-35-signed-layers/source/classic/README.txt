# Kissing-number improvements in dimensions 32–35

26 October 2026. Explicit configurations, checked by two separate exact programs.

| Dimension | Public comparison | This construction | Gain |
|---|---:|---:|---:|
| 32 | 352,320 | **352,800** | **480** |
| 33 | 370,240 | **370,368** | **128** |
| 34 | 390,980 | **390,982** | **2** |
| 35 | 415,308 | **415,310** | **2** |

The comparisons use Qihong Lin's October 2026 support codes, downloaded from
[Brouwer's live table](https://aeb.win.tue.nl/codes/Andw.html) on 26 October 2026:
1,713, 1,852, 2,013, and 2,202 supports, respectively. The standard construction
gives `131072 + 128*M + 2*n*(n-1)` points. These comparisons are newer than the
older values in the surrounding workspace. Source URLs and hashes are in
`source/baselines.json`.

Run the self-contained verifier, requiring only Python 3.9+:

```sh
python verify.py
```

It reconstructs the dense code, verifies its distance, checks the signed
replacements, and covers every pair of geometric layers by exact integer
inequalities. The four certificates verified in about three seconds on the
development machine. `verification.json` records the results and certificate
hashes.

A separately written checker expands the middle layers, roots, dense code,
and added heads into integer coordinate arrays. It requires NumPy:

```sh
python verify_coordinates.py
```

Its output is `coordinate_verification.json`. Neither verifier runs a search
or trusts the recorded verification output. `proof.md` gives the complete
construction and explains the normalization, including the two different
vector lengths used in dimensions 34 and 35.

Dimensions 32 and 33 use signed replacements of a few complete octad layers.
Dimensions 34 and 35 instead use a distance-six coset of the dense code and
replace two roots by four explicit heads. All four are lower bounds; no
optimality or exhaustive literature certification is asserted. The direct
fetch of Cohn's live table failed during this run; the comparison above is
against the successfully downloaded public code witnesses and the newer
benchmarks they imply.

Files:

- `certificate32.json` through `certificate35.json`: complete finite inputs.
- `verify.py`: standard-library exact checker.
- `verify_coordinates.py`: independent integer-coordinate checker.
- `proof.md`: formulas and pairwise verification argument.
- `verification.json`, `coordinate_verification.json`: recorded successful checks.
- `source/`: original public inputs, provenance, and the signed-support search
  and exporter. See `source/README.md` for reproduction commands.
- `SHA256SUMS`: hashes of all packaged files except the manifest itself.

The dense generator is the Cheng–Sloane `[32,17,8]` code, taken from the
already archived Qiushi certificate (`top_code.json`, revision
`bbb7aafd78e7c04d7e88260998355b5f3fc40d5c` of
[Qiushi-Engine-Kissing-Number-Research](https://github.com/Oxelra-AI/Qiushi-Engine-Kissing-Number-Research)).
It is reverified here. Layered binary constructions are due to
[Edel, Rains and Sloane](https://arxiv.org/abs/math/0207291);
signed replacements also appear in
[Yang et al.](https://arxiv.org/html/2609.35051v1).
The new claims concern these larger explicit configurations.
