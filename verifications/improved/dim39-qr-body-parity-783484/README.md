# 39D: 783,484 kissing points

An exact **+16** improvement over the [preceding 783,468-point construction](../dim39-qr-three-coordinate-heads/), not a proof of optimality.

The same 84 affine QR fibres and three physical head coordinates are retained. A simultaneous embedding changes 22 of the 36 physical body slots; as a permutation of the dense body coordinates (together with an allowed fixed body sign change), this preserves every internal dense angle. Short dual dependencies then align some sparse support projections with the body fibres. The five octad supports indexed in the certificate use odd sign parity; the other supports use even sign parity. Indices are zero-based in the support list **after** the seven inherited replacements. All incompatible octad signs are removed by the exact projected affine-fibre checks.

| Layer | Points |
| --- | ---: |
| Dense QR fibres | 344,064 |
| Signed octads | 436,448 |
| Roots | 2,908 |
| Stars | 64 |
| **Total** | **783,484** |

Run `python verify.py` in this directory. The verifier uses Python integers: it checks the shortened QR kernel, all 3,486 fibre-pair Hamming minima and exact angles, all retained sparse rows, the distinct-support intersection bound, and the root/star conditions. It reproduces 1,824 deleted octad signings, down from 1,840 in the preceding certificate. The three copied input files are the same Git blobs as in the [pinned predecessor](https://github.com/alexlegeartis/KissingNumbers/tree/d21336d488b917059c216632f10f6630d434a65c/verifications/improved/dim39-qr-three-coordinate-heads).

The QR source is [Yang et al., Appendix B](https://arxiv.org/html/2609.35051v1). The original 3,424 octad supports are from [Qihong Lin](https://arxiv.org/abs/2611.03682), CC BY 4.0, with the seven inherited support replacements in `base_certificate.json`.
