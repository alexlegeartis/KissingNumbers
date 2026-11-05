# Galois-ring trace codes and nonlinear translate search

This is a constructive algebraic backend, not a new kissing-number record. It builds Teichmuller representatives in GR(4,m), computes the Frobenius trace, synthesizes a Z4 generator, and applies the binary Gray map. The Gray image is nonlinear over F2: neither the search nor verification assumes that binary XOR is ring addition.

`galois.kerdock` independently reconstructs and exhaustively checks the Z4 Lee-weight spectrum. For m = 3 it generates the classical nonlinear (16,256,6) Nordstrom–Robinson code; for m = 5 it generates a nonlinear (64,4096,28) Kerdock code. The Lee-to-Gray Hamming isometry makes the enumerated minimum Lee weight an exact all-pair binary minimum-distance certificate.

`galois.find_translates` searches XOR-translates using the exact binary difference set and verifies every selected cross-distance. On the 16-bit control it finds eight layers, giving 2,048 binary sign points; the optional D16 root layer adds 480 and gives 2,528 points in dimension 16. That is below the known 4,320.

`galois.evaluation_search` searches additional order-two or arbitrary non-Teichmuller evaluation points. In one reproducible 96-coordinate example it finds a (96,4096,38) binary trace code. Its full spherical embedding has 4,096 points and exact rank 96, with all 8,386,560 pair inequalities independently checked. The repository's 96D record is vastly larger, so this is an algebraic capability and a search seed, not an improvement. The 16 extra evaluation points are chosen from the finite set of non-Teichmuller ring elements; no optimality is claimed.

`galois.translate_search` bridges the resulting Z4 generator to the general binary-translate proposal engine. Heuristic cross-checks are proposals; `galois.find_translates` returns an independently replayable exact certificate only after exhaustive verification. This backend does not reproduce the far larger Duursma–Greferath–Litsyn–Schmidt (96, 2^37, 24) code and is not a general discovery algorithm for arbitrary nonlinear codes.

Run from the repository root:

    python research/galois_kerdock_20261105/benchmark.py

`galois.random_translates` extends the same idea to binary lengths up to 128, without enumerating the whole ambient cube. It samples binary shifts, computes the **complete** XOR difference set of the nonlinear parent, and checks every accepted pair of translates exactly. The optional `wide96_benchmark.py` uses the searched (96,4096,38) seed and 100 discovered shifts, plus D96 roots. Its independently replayable result has 427,840 points in dimension 96, well below the repository's existing 142,162,362,560-point construction. This benchmark requires larger explicit pair, matrix and work budgets and is not a new record.
