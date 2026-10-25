# Support changes and root replacements in dimension 39

Checkpoint: **25 October 2026**.

| Dimension | Previous benchmark | This construction | Gain |
|---|---:|---:|---:|
| 39 | 768,916 | **768,924** | 8 |

The previous method is Qihong Lin's [ERS layers with an enlarged support code](https://github.com/QH-Lin/kissing-numbers-33-39/tree/0c4e2bd486061b3fd55e128ac3497bfb12b73f07).

Four cardinality-preserving changes to the weight-eight support family permit 64 new vectors of shape (4,2,2,2,2), replacing 56 roots. The verifier enumerates the dense binary code differences, checks all support intersections and every affected pair type with Python integers. The dense and middle layers retain 327,680 and 438,272 points.

```bash
python final.py
python verify.py
```

The wrapper checks the complete SHA-256 manifest and replays the standalone
proof in a temporary directory. The proof uses only Python's standard library.
`source/proof.txt` gives every angle argument and distinctness proof;
`source/certificate.json` gives the new points and exact transformation data.
`source/verify.py` is independent of search code, and `source/verification.json`
preserves the original receipt. `source/source/` includes all parent data and
`THIRD_PARTY.txt` attribution and licenses. `final.py` derives the count.
`data/SHA256SUMS.txt` covers every package file. No optimizer, search cache,
download or numerical tolerance is needed to verify this lower bound.
