# Independent norm-four layers in 82D and 83D

| Dimension | Verified lower bound | Previous repository bound | Increase |
|---:|---:|---:|---:|
| 82 | 6,269,102,538 | 6,269,087,562 | 14,976 |
| 83 | 6,278,834,796 | 6,278,767,980 | 66,816 |

We retain the entire [signed-apron construction](../dim82-83-gamma72-signed-apron/),
including the 6,218,175,600 contacts of [Nebe's Γ72](https://arxiv.org/abs/1008.2862)
and the original [Ganzhinov](https://arxiv.org/abs/2207.08266) and
[Bianchi et al.](https://arxiv.org/abs/2606.10402) tail coordinates.

**82D.** The old mixed tails have direction coordinates in a basis
weighted by 1 and 3. The new directions use the reciprocal ordering
of the weights (3 on the first five coordinates and 1 on the last
five). First take the 40 signed D5 roots on the first five coordinates.
A round-robin edge partition splits them into five sets of mutually
nonacute directions. The next 64 directions are 32 antipodal pairs
from the integer weighted norm-60 shell, encoded in witness.json.
All directions are separately normalized to squared norm four.

**83D.** Start with the 208 signed D11 roots compatible with the older
tails. Divide ten of the coordinates into five pairs. A tail has
entries ±1 on four whole pairs and zero on the fifth pair and the
unused coordinate. For a fixed omitted pair, there are four relative
sign bits. View them as an edge of the five-cube parallel to the
omitted coordinate. The 16 chosen edges form a perfect matching, so
different omitted pairs cannot agree on all common parity bits.
Parallel selected edges differ in at least two of their four fixed
bits, ruling out a single-sign discrepancy. Each edge gives 16 tails,
split into two sets of eight by parity of the first signs. Tails in
either set have nonpositive dot product. This gives 256 new octad
directions and 464 directions in all, after normalization.

Assign a distinct nonzero B-only phase to each set of nonacute tails.
Each phase has 144 horizontal norm-four vectors, in three separate
Leech-block orthogonal frames. The squared separation of two distinct
phases is at least four; it is at least six from any old mixed phase.
The new tails have dot at most two with old mixed tails and with all
new tails in other phases. Hence cross distances are at least eight.
At equal phase, distinct horizontal vectors have squared separation
at least eight, and distinct tails also do by nonacuteness. Old pure
tails are twice their old mixed tails, so their dot bound is four;
the Γ72 equator clears by the B-only phase depth. The old proof
checks the exact underlying Leech parity partition and coset minima.

Run `OPENBLAS_NUM_THREADS=1 python verify.py` with NumPy and SymPy and
assertions enabled. The pinned predecessor directories must be present
as siblings, or pass their parent directory as the first argument.
The proof uses exact integer inequalities, including exact sign tests
in the quadratic tail fields.
