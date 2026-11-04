# Validation of the initial implementation

The machine-readable receipt is [validation.json](../benchmarks/validation.json).
All results here are initial controls or bounded searches, not new SOTA claims.

- Full package suite: **55 tests passed**, 22.861 seconds.
- Optimized-Python core/coset/negative-replay suite: **25 tests passed**, 0.330 seconds. The DGLS
  assertion-based legacy proof adapter explicitly refuses optimized mode.
- Exact end-to-end controls: 510 points in 10D, 1932 in 14D, 352800 in 32D,
  and the implicit DGLS 137438953472-point sign layer in 96D.
- Numerical generator-based discovery: 510 in 10D; 6720 in 28D. Neither is a new
  record, and neither numerical receipt is promoted to an exact certificate.
- Independent extension replay: every one of 262144 linear extensions for each of
  two fixed Z4 modules excluded by 11/10 low-weight witness constraints. This
  statement does not apply to other modules or constructions.

Measured end-to-end times (including final proof replay) were 1.006s, 1.838s, 0.636s, 38.152s
for the 10D, 14D, 32D and 96D exact controls. Cumulative process peak RSS was
77840KiB. These are executor observations, not promised runtimes.

[Exact control receipts](../benchmarks/receipts/) preserve programs, computed
facts, proof summaries and source hashes. Run `python -m discovery verify PATH`
to recheck one. Receipt hashes detect alterations, not authenticated provenance.
The compressed 32D control performs 1826076 support/sign checks; the 96D code
remains implicit and uses two independent lower-ring/coset proofs.

The repository-wide audit passed all checks. It is run separately with
`python audit.py --repository-only`. New experiments do not register any result,
modify a sealed proof package, or change the 70 existing registered values.
No commit, push, Actions run, software installation or publication was performed.

The independent [final review](FINAL_REVIEW.md) found and verified fixes for three defects.


Later integration evidence is recorded in `research/discovery_campaign_20261104/FINAL_STATUS.md`, `REVIEW.md`, and `discovery/benchmarks/planner_integration_20261104.json`. The earlier measurements above remain historical snapshots.
