# A 33D neutral signed-line exchange

This gives another exactly checked 370,368-point coordinate configuration;
it is **not a new lower-bound record** and is not claimed to be nonisometric
to the source. It replaces 67 of the 320 exceptional antipodal lines with 67
new lines while retaining the full octad fibers, dense sign code and roots.
The changed lines were found in a pointwise ruin-and-repair search over all
signs of 23 new octad supports that have no full-fiber blockers.

From the repository root, run:

    python research/33d_neutral_signed_exchange_20261105/verify.py

The delta is relative to `verifications/improved/dim32-35-signed-layers/source/classic/certificate33.json`.
The verifier reconstructs the complete compressed construction and runs the
native exact signed-support checker; no search result or numerical tolerance
is used as the proof. The certificate does not claim a gain or priority.
