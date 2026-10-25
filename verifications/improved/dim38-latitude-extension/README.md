# Rational latitude extension in dimension 38

Checkpoint: **25 October 2026**.

| Dimension | Previous benchmark | This construction | Gain |
|---|---:|---:|---:|
| 38 | 605,322 | **605,346** | 24 |

The previous method is Qihong Lin's [Leech caps with an equatorial extension](https://github.com/QH-Lin/kissing-numbers-33-39/tree/0c4e2bd486061b3fd55e128ac3497bfb12b73f07).

Six new head lines give 24 points at a new latitude. An exactly orthogonal rational rotation fixes these heads and retains all 1,932 axis points. The verifier rebuilds the Golay Leech shell, checks the 589,680 caps, all 13,710 equatorial points, and every new angle type. Axis/cap products use exact integer interval bounds with proved overflow limits.

```bash
python final.py
python verify.py
```

The wrapper checks the complete SHA-256 manifest and replays the standalone
proof in a temporary directory. NumPy is required.
`source/proof.txt` gives every angle argument and distinctness proof;
`source/certificate.json` gives the new points and exact transformation data.
`source/verify.py` is independent of search code, and `source/verification.json`
preserves the original receipt. `source/source/` includes all parent data and
`THIRD_PARTY.txt` attribution and licenses. `final.py` derives the count.
`data/SHA256SUMS.txt` covers every package file. No optimizer, search cache,
download or numerical tolerance is needed to verify this lower bound.
