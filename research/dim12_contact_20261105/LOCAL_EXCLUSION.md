# An exact local obstruction around a 12D 842-point near-miss

This certificate excludes a specific small neighborhood. It does **not** exclude
842-point kissing arrangements in dimension 12, prove global optimality of a
near-miss, or prove that its contact graph must change.

Let the supplied coordinate rows define mathematically normalized unit rays
`x_i`. Choose nonnegative rational weights `lambda_ij` on selected pairs,
summing to one. Define

\[
 g=\sum_{i<j}\lambda_{ij}x_i\cdot x_j,\quad
 a_i=\sum_j\lambda_{ij}x_j,\quad
 h_i=a_i\cdot x_i,\quad f_i=a_i-h_ix_i.
\]

Here all supported pairwise dots are positive, so every `h_i >= 0`, and
`sum_i h_i = 2g`. No equilibrium assumption `f_i=0` is made.

For any other unit vectors `y_i`, write `delta_i=y_i-x_i`. The unit-length
identity gives `x_i . delta_i = -||delta_i||²/2`. Therefore

\[
 \sum_{i<j}\lambda_{ij}y_i\cdot y_j
 =g+\sum_i f_i\cdot\delta_i
  -\tfrac12\sum_i h_i\|\delta_i\|^2
  +\sum_{i<j}\lambda_{ij}\delta_i\cdot\delta_j.
\]

If every `||delta_i|| <= r`, Cauchy–Schwarz yields

\[
 \max_{i<j}y_i\cdot y_j
 \ge g-r\sum_i\|f_i\|-r^2(1+g).
\]

A lower bound greater than `1/2` rules out a kissing arrangement in this ball.
This conclusion remains valid when the fitted forces are merely small; the
entire force residual is explicitly paid for in the bound.

## Exact arithmetic

`local_certificate.py` interprets every stored binary64 coordinate as its exact
dyadic rational. It encloses the normalized rays using integer square roots
and outward rational rounding with denominator `2^80`. Fitted numerical
weights are rounded to nonnegative integers at scale `2^48`, then divided by
their actual integer sum, giving exactly normalized rational weights.

All dot, force, and force-norm bounds use integers and rational arithmetic.
The force norm is bounded with an upward-rounded integer square root.
The final decision uses

\[
 g_{\rm lower}-rF_{\rm upper}-r^2(1+g_{\rm upper})>\tfrac12
\]

as an exact rational inequality. Numerical optimization finds useful weights,
but is not trusted to establish any inequality in the certificate.

The supplied files preserve the actual coordinate and weight inputs. They can
be checked independently of PyTorch, CUDA, the numerical stress fitter, or any
claim that its fit converged.

## Certified radii

For the literal rays in the stored reconstruction of the published near-miss,
`evidence/stress_published842.npz`, the certificate gives

* `r = 6121/250000 = 0.024484`;
* weighted average dot about `0.5008999065`;
* the exact lower bound exceeds `0.50000005728`.

Thus a valid 842-point arrangement must move at least one of the matched points
by more than this chord distance (about 1.4 degrees). The statement also applies
after any matching and rotation that places all the candidate points within the
specified distances: that would itself supply the forbidden labelled `y_i`.

For the improved near-miss, `evidence/stress_improved842.npz`, the independently
checked radius is `24409/1000000 = 0.024409`, with a lower bound exceeding
`0.50000004163`. This smaller radius reflects the particular fitted weights and
residuals; it is not a claim about which configuration has the larger true
feasibility neighborhood.

The rational values and input hashes are in `published_local_certificate.json`
and `local_certificate.json`. Replay with:

```sh
taskset -c 0-126 python - <<'PY'
from local_certificate import certify
from pathlib import Path
for label in ['published', 'improved']:
    result = certify(Path('evidence') / f'stress_{label}842.npz',
                     Path(f'rechecked_{label}.json'))
    print(label, result['radius'], result['diagnostic_lower_bound'])
PY
```
