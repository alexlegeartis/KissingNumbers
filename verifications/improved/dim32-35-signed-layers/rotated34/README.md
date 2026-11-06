# 390,992 kissing points in dimension 34

This 6 November 2026 construction improves the preceding 390,990-point repository certificate by two points. It retains the 2,013 public weight-eight supports, changes the dense coordinate placement, switches one star parity and rotates a heavy head rationally. It is an explicit finite construction; no geometric minimum, population or numerical optimization premise is inherited.

## Construction

The layers contain 131,072 dense points, 257,664 middle points, 64 stars, four heavy heads and 2,188 roots, totaling **390,992**.

Place the dense [32,17,8] code and its recorded distance-six coset on the coordinates other than 0 and 2. Keep all 128 even sign patterns on each recorded eight-support. Use the four preceding star families, but make the leaf parity odd for the family centered at 0 with leaves 2,4,14,29; the other three families retain even parity. Delete exactly the preceding 28 root lines. In particular, the opposite-sign (0,2) root line is already deleted, while the formerly deleted (32,33) line is restored.

Let h1 and h2 have entry 1 away from coordinates 0,2, and tail entries

- h1: (123/85, -339/85),
- h2: (-3, 3).

Both have squared norm 50. Include ±(4/5)h1 and ±(4/5)h2. The rotation makes h1 compatible with the switched star parity, so both heavy lines share an already necessary root deletion. This restores two roots without losing another point.

All coordinates are rational. Multiplication by 425 gives integer coordinates of squared norm **5,780,000**, with every distinct inner product at most **2,890,000**.

## Exact verification and reconstruction

The parent package's `verify.py` replays both routes and rejection controls in temporary copies:

- `rotated34/verify.py` enumerates the dense code, checks all 2,025,078 support pairs, and proves the interaction bounds using exact integer parity and coordinatewise calculations.
- `rotated34/verify_coordinates.py` independently expands all 390,992 vectors, checks all 15 interaction types, verifies uniqueness and exact rank 34, and imports neither the symbolic verifier nor historical verification code.
- `rotated34/negative_controls.py` requires both routes to reject duplicate supports, a wrong head norm, an undeleted incompatible root, a false count and a duplicate star. `negative_control_results.json` records ten rejections.

To reconstruct the complete common-norm integer witness outside the package:

```bash
python rotated34/verify_coordinates.py --coordinates /tmp/kissing34-390992.npz
```

The archive contains a `coordinates` array of shape (390992,34). A large binary witness need not be downloaded: `certificate.json` and the verifier reconstruct it without search. `verification.json` and `coordinate_verification.json` record the two independent replays. The certificate SHA256 is `36a9017be6de898d16bdfdad1b686db1d991a9be8a77a804134d770794a5f7ea`.

## Comparison and scope

The public Lin support-code/ERS construction supplies 390,980 points; the preceding repository certificate has 390,990. Its source code and root-star construction remain unchanged in the parent package. This result adds two to that repository record, twelve to the pinned public witness. The existing published comparison is retained; this computation is not a new comprehensive worldwide literature audit. The manuscript is intentionally unchanged.
