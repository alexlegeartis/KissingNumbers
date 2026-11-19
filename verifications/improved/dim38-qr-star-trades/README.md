# 38D: 652,044 by coordinated star trades

Start with the [652,028-point QR construction](../dim38-qr-free-offsets/).
Keep its 253,952 dense points and 395,264 signed octads. Add 144 stars
and remove exactly the 128 signed roots incompatible with them,
leaving 2,684 roots and **652,044 points**.

A star has one coordinate equal to ±2 and four coordinates equal to ±1,
so its squared norm is 8. The four columns of `stars.txt` are the
zero-based center, the leaf mask, the center sign, and the negative-leaf
mask. A star alone conflicts with just four signed roots, one for each
center-leaf pair with matching signs. The gain comes from choosing stars
whose forbidden roots overlap. Selecting them also requires pairwise
star products at most 4.

For each star and retained octad support, the sum of the magnitudes on
the support is at most 4, giving raw product at most 8 with the signed
octad at squared norm 32. The QR dense points are handled by the
coordinatewise magnitude bound. The standard-library replay checks
both and exactly counts the union of the forbidden roots; the NumPy
audit explicitly tests every star against all 2,812 original signed
roots and enumerates every star-star pair independently.

```sh
python verify.py
python audit.py
```
