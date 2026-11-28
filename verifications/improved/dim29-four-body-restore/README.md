# 232,244 kissing points in dimension 29

This restores a missing antipodal tail pair, adding **two** points to the [232,242-point shared-tail splice](../dim29-shared-tail-reflection/). The result is **14 above [Rustem Takhanov's 232,230-point v0.2 witness](https://github.com/k-nic/KissingNumbers/tree/77f0860bf93a745140e6e25481d3f60773b11441/version2)**. The construction extends Takhanov's Leech-midpoint scaffold; it retains his complete 164,362-point equator and all auxiliary points.

## Four-body colour switch

Write a paired midpoint fibre as \((b,\pm t)\), with \(b\) the 24-coordinate body and \(t\) the 5-coordinate tail. For two orthogonal tail lines of squared norm 32, equally coloured bodies require body dot at most 32; differently coloured bodies can have body dot up to 64. Above 32, colours in a connected component have to switch together.

Take the lines represented by \(B=(4,0,0,0,-4)\) and \(C=(0,4,0,4,0)\). In the pinned predecessor, the missing target body with proposed tail \(\pm B\) is blocked by one paired body of colour \(B\). The forced two-colour component contains four paired bodies, two of each colour. Swapping both colours on the whole component clears the target and restores its two rays. [The literal eight removals and ten additions](edits.json) are enough to verify the construction without relying on the search.

## Exact replay

First generate the predecessor file `coordinates29_232242.npz` using the [two previous replays](../dim29-shared-tail-reflection/#exact-replay), beginning with Takhanov's pinned original coordinates. Then, from the repository root:

```sh
python verifications/improved/dim29-four-body-restore/verify.py \
    coordinates29_232242.npz --output coordinates29_232244.npz
```

The verifier checks the predecessor's canonical hash before the edits. All 232,234 unchanged rows inherit their pair constraints. The ten changed/new rows have squared norm 128; every changed–retained and changed–changed inner product is checked in signed integer arithmetic, with both maxima exactly 64. The canonical final SHA-256 (lexicographically sorted little-endian int16 rows) is `535725b88ca163f30194802929e58046e940ad041932341abf0beab5ce91cffa`.
