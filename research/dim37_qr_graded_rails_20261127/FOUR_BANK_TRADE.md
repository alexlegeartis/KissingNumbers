# 37D QR: exact whole-cubic tie via four deep-hole banks

This small independent certificate concerns the same `dense_parent38.json` (see the companion graded quotient obstruction). At the two deepest ordinary QR cosets relative to +4 dense rows, the complete four distance tuples in grade order (+4,+1,−4,−1) are `(10,1,9,10)` and `(10,9,1,10)`. Their sign-reflected graded cosets use the same XOR grade-2 offset. The four new bank generators are:
```
12518842  +4
 8379313  -4
14022956  +4
 1386279  -4
```
Syndrome free positions are ordinary indices 0…23 in this pinned elimination, so the decimal numbers also give literal full ordinary mask representatives with indices 24…35 zero. Each bank contains all 4,096 masks given by this representative plus the projected zero-grade QR kernel, putting amplitudes ±2 on all 36 ordinary sites and its listed scalar head in coordinate 36. Its squared norm is 160.

Compare against every old bank (11 mother centers, 4 graded parts each) by exact min ordinary Hamming distance among the 4,096 binary kernel words. Each of the four new banks conflicts with exactly one old bank at distance 1: respectively tail grades 1,3,2,0, **all in zero-based mother center 7**. Removing that old mother center releases all of them. The four new banks have mutual minimum distances 10 or 11 as appropriate, and each independently meets every vector in every remaining old dense bank. Within one new bank the minimum Hamming distance is at least 12 (zero-grade kernel); the norm-160 gate has slack.

The changed banks also meet every nondense point, without assumptions on their phases: their raw dot with any norm-32 octad is at most 32; with an ordinary root at most 16; with a radical six-support at most 24+8√2<√(160·32)/2; and with either unit marked pole at most 4<√160/2. This is an exact 16,384-for-16,384 exchange, remaining at the same existing 544,666 bound, not an improvement. It exposes a specific whole-center redesign rather than a template parameter sweep.
