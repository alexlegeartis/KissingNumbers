# A 43D one-for-two head exchange

This package increases the repository's certified 43D construction from
2,919,782 to **2,919,783** points. It removes one point in the explicit tiny
layer and inserts two rational-head points on the same tiny body. All 118
complete families and five pure heads are retained. The proposal came from a
finite search coupling existing head directions with the 30,240 tiny bodies;
it was not discovered by the unified engine.

From the repository root:

    OPENBLAS_NUM_THREADS=1 python research/43d_head_exchange_20261104/replay.py

The script checks the baseline archive hash, constructs the delta in a
private temporary copy, and runs both existing exact verifiers: one in the
ambient 48-coordinate model and one in the independent rank-40 coefficient
model. The second replay also checks that the candidate is exactly the stated
one-for-two exchange. Both use integer comparisons for every angle decision.
The complete-family population identities and the P48n minimum/design inputs
are inherited from the existing package; this replay does not freshly
enumerate all large shells or independently prove those inputs.

The `delta.json` indices are zero-based, relative to the original
`partial43/explicit_points_43.npz` at the pinned commit. A numerical search
or optimizer is not needed to reproduce or validate the certificate. The
current literature priority should be checked separately; the comparison here
is against the repository's 2,919,782-point certificate.

## Discovery replay

    OPENBLAS_NUM_THREADS=1 python research/43d_head_exchange_20261104/search.py --output /tmp/k43-delta.json

The search crosses the 30,240 tiny bodies with the 305 distinct rational
heads in the incumbent. It screens against the complete families and the
incumbent tiny layer, groups candidates by their unique deleted incumbent,
and looks for two mutually compatible replacements. Integer arithmetic
checks the proposed pair against all retained tiny points, the complete
families and the pure heads before emitting it. A floating-point screen is
used only to generate candidates and may miss feasible exchanges; it is not
a certificate or an exhaustive optimality claim. On the development run it
reproduced `delta.json`; the independent full checks remain `replay.py`.
