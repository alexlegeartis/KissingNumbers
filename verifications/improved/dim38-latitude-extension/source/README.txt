VERIFIED LOWER BOUND: tau_38 >= 605346
25 October 2026

Run:
    python verify.py

Requires Python and NumPy. Verification takes a few seconds and uses
integer arithmetic and exact rational interval bounds. The program
rebuilds and checks the existing cap/equator configuration, the rational
axis rotation, and every pair type involving the 24 new points.

Read proof.txt for the construction. certificate.json and source/ contain
all required data. verification.json records the successful exact check.

This improves the verified public benchmark 605322 by 24 points. It is
a local exact certificate, ready for independent review; no external
submission or record-table acceptance is claimed.
