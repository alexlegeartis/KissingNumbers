VERIFIED LOWER BOUND: tau_39 >= 768924
25 October 2026

Run:
    python verify.py

Only the Python standard library is required. The exact integer verifier
takes about one second, including the dense-code and support checks.

Four neutral support changes allow 64 new vectors to replace 56 roots.
Read proof.txt for the construction. certificate.json and source/ contain
all required data. verification.json records the successful check.

This improves the verified public benchmark 768916 by eight points. It is
a local exact certificate, ready for independent review; no external
submission or record-table acceptance is claimed.
