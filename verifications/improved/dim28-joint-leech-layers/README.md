# Dimension 28: jointly selected Leech blocks and auxiliary layers

Checkpoint: 24 October 2026. **205,232 points**, improving Kevin Li's
205,172-point configuration by **60**.

The construction retains Li's 512 odd Leech heads in each of eight blocks
and jointly selects compatible even heads and an auxiliary layer over
dual-D4 directions. Four blocks have 524 heads and four have 528, for
4,208 owners. The count is

    192,352 equator + 3 * 4,208 caps + 256 auxiliary = 205,232.

Returning 16 owners to the equator and removing their caps costs 32 net
points; the auxiliary layer gains 92. The earlier 204,896-point frame-layer
construction remains in [`../dim28-norm8-frame-layer/`](../dim28-norm8-frame-layer/).

## Sources and credit

The lifting method, odd blocks and starting configuration are from
[Kevin Li's paper and certificates](https://github.com/Felpix-Studios/kissing-number/tree/b823640f4e98134b9d55b27288132ee954a8affa/preprints/leech-lift-d28-d31).
The joint selection and new certificate were obtained in this project.
The source certificate hash and input archive hash are recorded in
`source/README.txt` and `source/D28-205232.provenance.json`.
The unmodified all-pairs checker and its Apache 2.0 license are retained.

## Reproduce

```bash
python final.py
python verify.py
```

Python 3.9+ and NumPy suffice. The wrapper checks all file hashes and
replays the independent Python checker in a temporary copy. It rebuilds
the Golay Leech shell and checks its complete partition, every head-set
angle, all cap-tail inequalities and every auxiliary interaction with
unbounded integers in Q(sqrt(2),sqrt(3)). Equator/equator pairs use the
standard kissing property of the generated Leech minimal shell. No
floating-point tolerance enters acceptance.

The supplied all-pairs log records **21,059,984,296 pairs**, zero violations,
and PASS with Kevin Li's unmodified C++ checker. To repeat that independent
check, compile `source/check_kissing.cpp` as described in `source/README.txt`
and pass `source/D28-205232.txt`. The structured replay is the default job.

The search fixes the odd blocks and uses a finite candidate pool. Its
feasible objective is 480 and its solver upper bound is 530 within that
model; the construction family remains open.

## Files

- `final.py` derives the count from the certificate's class sizes.
- `verify.py` checks integrity and runs the isolated exact replay.
- `data/SHA256SUMS.txt` covers this package except itself.
- `source/D28-205232.txt` is the explicit coordinate certificate.
- `source/D28-205232.provenance.json` records the selected heads and inputs.
- `source/verify.py` is the independent exact Python checker.
- `source/check_kissing.cpp` is Li's unmodified all-pairs checker.
- `source/THIRD_PARTY_LICENSE.txt` preserves the Apache 2.0 license.
- `source/README.txt` and `source/SHA256SUMS.txt` preserve the standalone
  release's instructions and byte inventory.
- `source/D28-205232.verification.json`, `source/structured28_205232.log`
  and `source/all_pairs28_205232.log` preserve its recorded verification runs.
