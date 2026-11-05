# 13D exchange over Q(sqrt(3))

This finite-pool experiment starts from a supplied 1,154-point 13D code with
coordinates a + b sqrt(3) and squared norm 16. Its 1,728 candidate directions
come from a separate sparse algebraic shell and each conflicts with at most six
source points. Both the source and pool are included in `pool.npz`; this
package does not reconstruct their discovery. The target is a gain of at least
one point. The input is replayed exactly, and any solver proposal is checked
independently using integer arithmetic in Q(sqrt(3)).

    OPENBLAS_NUM_THREADS=1 python research/quadratic_exchange13_20261105/run.py --seconds 120

No improvement is asserted. A timeout or failure to find a proposal is not a
proof of impossibility. A candidate above 1,154 still requires independent
verification and comparison with the then-current literature before a record
claim. The solver explores this fixed finite pool, not all 13D configurations.
