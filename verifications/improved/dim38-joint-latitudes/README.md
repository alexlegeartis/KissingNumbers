# Joint latitude constructions in dimension 38

Checkpoint: **6 November 2026**. Exact lower bound: **611,638**, an improvement of **116** over the preceding repository certificate **611,522**.

589,548 retained caps + 13,710 equator points + 1,568 rotated flat axes + 28 coordinate axes + 6,272 Fano-grid points + 512 orthogonal Hamming points = **611,638**.

The [Hamming construction and exact proof](hamming/README.md) replace the 384 cube-edge points by 512 points on full seven-coordinate heads. Eight orthogonal body lines are colored by the binary Hamming syndrome, with a new body/tail norm split. Four old cap lines are restored and six new ones are removed: 12 additional cap deletions give a net gain of 116.

```bash
python final.py
python verify.py
python hamming/construct.py --check
python hamming/negative_controls.py
```

`final.py` derives the count. `verify.py` checks the complete `data/SHA256SUMS.txt` manifest and replays the unchanged historical source, both current exact verifiers and four corruption controls in a temporary directory. `hamming/verify.py` checks all 130,816 new pairs and every cross-layer condition. `hamming/verify_independent.py` independently rebuilds the source data and checks cap classes, equator pairs, all retained-grid pairs, restored-cap interactions and the new layer. Its retained cap-cap argument uses the source's elementary Golay/Leech minimum-norm proof, not an exhaustive Leech-pair enumeration or a new formal proof of that theorem.

Python and NumPy suffice. All geometric checks use integer products and sign-preserving symbolic squaring, without floating-point tolerances. The certificate, integer point data and receipts are under `hamming/`; its README lists each file.

## Preserved history and attribution

The [28 October Fano-grid and cube-edge construction](refill/README.md) reached 611,522 from 610,394, adding 1,248 latitude points and deleting 120 caps. Its 6,656-point population saturates that particular 63-support model. The new full seven-coordinate supports and changed norm split lie outside that ceiling. `refill/verify.py`, `refill/verify_coordinates.py`, `refill/construct.py --check` and `refill/negative_controls.py` remain unchanged. The expanded historical replay checks all 22,147,840 latitude pairs.

The retained 26 October construction has 1,352 four-point latitude families and 610,394 points; its proof is `source/proof.txt` and its replay is `source/verify.py`. The earlier local rational construction has 605,346 points. [Qihong Lin's public Leech cap extension](https://github.com/QH-Lin/kissing-numbers-33-39/tree/0c4e2bd486061b3fd55e128ac3497bfb12b73f07) has 605,322 points. Historical inputs and attribution are unchanged. No new comprehensive priority survey or global optimality claim is made. The manuscript remains deferred.
