# Two shifted polar additions in dimension 65

Keep the self-derived cubic/sparse equator and both axial poles, giving
the [prior verified 413,526,914 points](../dim62-67-rm2-selfsparse/). Define a linear subcode H of
RM(2,6), imposing six parity checks on the 15 quadratic coefficients
and leaving affine terms unconstrained. The syndrome bit pattern is the
field sum of `a*b*(a+b)` for every active quadratic coefficient indexed
by distinct powers of two a,b. The field basis uses `x^6+x+1`. The
linear subcode has 65,536 words, minimum distance at least 24; both
properties are checked exhaustively, so field irreducibility or an APN
theorem is not needed as a premise.

Take the two affine translates of H whose cubic label and quadratic
syndrome are (275722,44) and (686811,62); use the first quadratic
representative of each syndrome in mask order. Map their 131,072 body
words to ±1 signs and append the common coordinate +sqrt(17). For the
second layer, change both representatives by the first representative
of syndrome 1 and append -sqrt(17).

The exhaustive checks are:
* Body difference against the entire original dense equator has
  Hamming distance at least 14 under every old cubic label, assessed
  directly against all 2^22 quadratic words.
* Same-side differences have Hamming distance at least 22, including
  unique labels and all subcode differences.
* Opposite-side differences have distance at least 8.

Each new point has squared norm 81. Against the original norm-64
equator, its inner product is at most 36 — via distance 14 for dense
signs, and at most 32 by the sparse coordinate bound. The poles attain
at most 8 sqrt(17) < 36. Same- and opposite-side products are at most
37 and 31, respectively, both below half of 81. Thus after individual
normalization all constraints hold and the union gives **413,789,058**
points in dimension 65.

```sh
g++ -std=c++17 -O3 -march=native lift65.cpp -o lift65
./lift65 witness.txt
```
