Verified lower bound: tau_39 >= 769180

Run: python verify.py
Dependency: NumPy (the certificate was tested with NumPy 2.0.2).
Do not use python -O. No search software is required.

proof.txt                 Complete proof and construction definition.
orbit_certificate.json    Exact integer representatives and binary group.
verify.py                 Complete integer orbit verifier.
orbit_verification.json   Output of the successful verification.
base_certificate.json     Sparse-layer support replacements and heavy stars.
verify_base.py            Inherited verifier for the base sparse layer.
base_proof.txt            Earlier construction's proof; dense layer superseded.
source/                   Original source data and attribution.

Count: 1281 orbits of order 256, plus 441244 sparse points.
Gain: +248 over the preceding local bound 768932;
      +264 over Lin's public benchmark 768916.
Local exact verification does not assert external review or acceptance.
