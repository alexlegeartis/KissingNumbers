# 232,242 kissing points in dimension 29

This adds **six** points to the exact [232,236-point two-colour witness](../dim29-kempe-tail-switch/) and is **12 above [Rustem Takhanov's 232,230-point v0.2 witness](https://github.com/k-nic/KissingNumbers/tree/77f0860bf93a745140e6e25481d3f60773b11441/version2)**, the competing value in the repository's dated comparison. The construction extends Takhanov's Leech-midpoint scaffold; it does not attribute that scaffold to this repository. Coordinates have squared norm 128 and distinct dots at most 64.

## Share the cost of conflicting rows

Write a 29D layer point as \((b,t)\in\mathbb R^{24}\oplus\mathbb R^5\). Reflecting just the tail gives \((b,-t)\). The missing tail antipodes considered in this search are mutually compatible: their originals were all in a legal configuration, and tail reflection is one common isometry.

Build a bipartite graph between the proposed tail antipodes and the occupied baseline rows, joining pairs with forbidden raw dot \(>64\). Inserting a set \(S\) requires deleting only its distinct old neighbors \(N(S)\), so the net gain is \(|S|-|N(S)|\). A matching-based search found **ten new rows sharing just four conflicting old rows**, yielding 232,242. The literal edits are in [edits.json](edits.json). The proof does not require assuming the search optimal or trusting a matching implementation.

## Exact replay

From the root of this repository, start with [Takhanov's pinned coordinates](https://github.com/k-nic/KissingNumbers/blob/77f0860bf93a745140e6e25481d3f60773b11441/version2/certificates/data/d29_coordinates.npz) and rebuild the existing Kempe witness. Then apply the literal edits:

```sh
OPENBLAS_NUM_THREADS=1 python verifications/improved/dim29-kempe-tail-switch/verify.py \
    d29_coordinates.npz --output coordinates29_232236.npz
python verifications/improved/dim29-shared-tail-reflection/verify.py \
    coordinates29_232236.npz --output coordinates29_232242.npz
```

The predecessor can also be replayed using its `--full` option. The new verifier checks the predecessor's canonical coordinate hash, the ten squared norms, every new–retained and new–new dot in signed integer arithmetic, and the new canonical coordinate hash. Both dot maxima are 64. Retained–retained pairs are exactly inherited from the pinned predecessor.

Canonical SHA-256, with rows lexicographically sorted as little-endian int16:
`ac8fb373b0321a42ce71ffbeef0dde75b9e74a3e4c8ebb807d24ba20d910d848`.
