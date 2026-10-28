# Joint latitude constructions in dimension 38

Checkpoint: **28 October 2026**. Exact lower bound: **611,522**.
The Fano-grid and cube-edge construction improves the preceding repository
bound **610,394** by **1,128** points.

589,560 retained caps + 13,710 equatorial points + 1,568 rotated flat axes
+ 28 coordinate axes + 6,272 Fano-grid points + 384 cube-edge points
= **611,522**.

The new latitude layer adds 1,248 points while deleting 120 caps.
[The construction and proof](refill/README.md) give an explicit small
recipe and a ceiling for its particular latitude support model.

```bash
python final.py
python verify.py
python refill/negative_controls.py
```

The wrapper first replays the unchanged historical 610,394-point
certificate in `source/`, then checks the current `refill/certificate.json`
with `refill/verify.py` and `refill/verify_coordinates.py`. The latter
independently reconstructs all original layer facts and checks every
distinct pair among the 6,656 expanded latitude points. It also runs
`refill/construct.py --check` against `refill/recipe.json`.
`refill/negative_controls.py` tests rejection of four corruptions.

The current CI receipts are archived in `refill/verification.json` and
`refill/coordinate_verification.json`.

Python and NumPy suffice. All mathematical checks use integers. The
SHA-256 manifest covers every package file other than the manifest itself;
replays run in temporary directories.

The retained 26 October construction has 1,352 four-point latitude
families and 610,394 points; its detailed proof is `source/proof.txt`.
The earlier local rational latitude construction has 605,346 points.
[Qihong Lin's public Leech cap extension](https://github.com/QH-Lin/kissing-numbers-33-39/tree/0c4e2bd486061b3fd55e128ac3497bfb12b73f07)
has 605,322 points. These historical inputs and attribution remain
unchanged. Gains above compare with the preceding repository certificate;
no new comprehensive priority survey or optimality claim is made.
