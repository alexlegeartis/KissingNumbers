# A 232,157-point kissing configuration in dimension 29

The configuration replaces 11 points of the [232,154-point parent](../dim29-takhanov-auxiliary-exchange/)
with 14 new points. It improves that verified lower bound by three and
the [pinned Takhanov comparison](https://github.com/k-nic/KissingNumbers/tree/dadcb2940f7c3c3bd478b0543931d3eef9f4bd65/version1)
by four.

## Construction

Work with integer numerators of squared norm 128, the same scale as the
parent. The new vectors are supported on the 13-coordinate section shown
as `coordinates` in `trade.npz` (0-based coordinate indices). Enumerate
all even integer vectors of squared norm 128 in this section and price
each ray by the set of parent points with dot product greater than 64.
For the blocked-point calculations, identical projections onto the
section were collected with their exact multiplicities. This is valid
because a ray in the section has the same product against a parent point
as against its projection. The discovery search found three compatible
trades and combined them, taking the union of the blocked indices only
once. The statement proved is the existence of the explicit final code,
not saturation of this or any larger search space.

The compact certificate records each primitive signed ray, its squared
norm, and the exact 11 deleted parent indices. The verifier scales all
rays back to squared norm 128, checks each against every one of the
232,154 source points in signed integer arithmetic, and checks all 91
pairs of distinct new points. Nondeleted parent pairs are inherited from
the pinned parent certificate. An independent all-pairs replay verifies
all 26,948,320,246 distinct unordered final pairs with binary32 matrix
products that are exact for these small integer inputs.

## Replay

From the root of the repository, using Python and NumPy:

```sh
python verifications/improved/dim29-even-grid-exchange/verify.py
python verifications/improved/dim29-even-grid-exchange/verify.py --full
python verifications/improved/dim29-even-grid-exchange/verify.py --output /tmp/coordinates29_232157.npz
```

The first command checks the exchange using the source bytes pinned in
`manifest.json`. The second checks all pairs without relying on the
parent's geometry proof. The third exports the complete coordinate code.
The full replay result is in `verification.json`.
