# Current tiny-point completion in43 dimensions

The4 November2026 certificate gives **tau(43)>=2,919,785**. It improves the
preceding2,919,782 by three, through one-for-two and two two-for-three exchanges.
All118 complete families and five pure heads remain unchanged. There are3450
explicit tiny points.

The generic engine generated the same final point set as the separate43D
search, with the last two steps in a different order. No delta or target count
was supplied to the operator. Original discovery and seeded engine reproduction
are separately recorded in `UPDATE_PROVENANCE.json` and research packages.

The current `delta.npz` is composed against `previous_points.npz` (3408 points).
`reproduce.py` checks exact equality with the registered arrays. The original
3447-point certificate and its delta, reports and notes are retained under
`history_2919782/`; its delta can still be checked using explicit `--points`
and `--delta` arguments to `verify_coordinates.py`.

Run from the repository root:

```sh
python verifications/improved/dim41-46-p48-heads/continuation43/reproduce.py
python verifications/improved/dim41-46-p48-heads/verify.py
python verifications/improved/dim41-46-p48-heads/continuation43/negative_controls.py
```

Both independent geometry replays use exact integer decisions. Complete-family
populations retain the inherited identities and P48n minimum/design inputs in
`../POPULATIONS.txt`. This update does not freshly enumerate the large shells
or assert optimality or worldwide priority.

The archived28 October notes describe the earlier continuous-head search.
Its optimizer, pool and search receipt remain historical inputs; the current
coordinate proof needs only the predecessor and the composed exact delta.

Historical optional search tools are `search.py` and `select_candidates.py`
(`continuation43/search.py` and `continuation43/select_candidates.py` from
the package root). Their pool and search receipt describe the28October
checkpoint; they do not reconstruct the current combined delta.
