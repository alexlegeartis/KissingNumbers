# Verified 35-dimensional kissing configuration: 415,314 points

This certificate improves the earlier local result of 415,310 by four points.
It improves the 415,308-point construction obtained from Qihong Lin's
2,202-support code, downloaded from Brouwer's table on 2026-10-26, by six.
This is a small improvement, not the substantial increase sought in this run.
The 34-dimensional result remains 390,982.

The certificate specifies all points exactly. Normalize each base-layer vector
by the square root of 32, and each added head by the square root of 51.
Both signs of every listed antipodal line are included.

Run from this directory:

```text
python verify.py
python verify_coordinates.py
```

The first checker uses only the Python standard library. The second requires
NumPy and independently expands the coordinate layers, checking 4,983,624
head/base dot products with integer arithmetic. Neither checker imports the
constructor or heuristic search code. Verification uses assertions; do not run
Python with `-O`.

To regenerate the certificate deterministically:

```text
python source/multi_head.py --baseline source/original_certificate35.json --output certificate35.json
```

`proof.md` explains the construction. `SHA256SUMS` identifies the packaged
files. The source input preserves the previous certificate, including the
source URL and SHA-256 digest of the downloaded public code. The corresponding
original public code file is included in `source/`.

This is a verified lower bound, not an optimality proof or an exhaustive
certification against every unpublished or independently published result.
No public registry or manuscript was changed.
