# Delsarte–Goethals spike layers in 60D and 61D

| Dimension | Predecessor | Net added contacts | New lower bound |
|---:|---:|---:|---:|
| 60 | 68,721,704 | 982,920 | **69,704,624** |
| 61 | 118,850,346 | 1,474,440 | **120,324,786** |

The starting configurations are the
[verified transferred support banks](../dim60-61-rm80-guard-transfer/).
The new construction keeps all their dense, affine-plane and guarded
sparse banks. Remove only the 120 ordinary coordinate poles, which
would otherwise conflict with the larger spikes. In 61D both tail
poles remain.

Assign a different published mother cubic translate to each of the
60 ordinary coordinates, using the first 60 in the pinned order.
Within the quadratic coordinate of each assigned translate, use the
1,024 quadratic forms of the explicit binary Delsarte–Goethals code
in dg.py. Its distance after arbitrary affine translations is at
least 24. The four marked sites form an affine plane, so the affine choices
project to the eight even-parity labels, each with multiplicity 16.
The same quadratic family has exactly 512 forms of each restricted
four-bit parity. Adding the fixed cubic merely exchanges the parity
classes if its restriction is odd. Each proposed label therefore
admits 512 quadratic forms and 16 affine choices per form, giving
8,192 contacts per label per ordinary spike position. The labels in the new layer are
absent from the old dense carrier, so a donor can never share its
entire unpunctured binary message with an old dense vector.

* **60D:** use omitted labels 2 and 3 on marks (0,20,41,61), spike
  amplitude 11/2, all other ordinary magnitudes one. The donor norm
  squared is 357/4. The old dense labels 0 and 1 may be at Hamming
  distance two from a donor label; the resulting exact angle check
  passes even if the common sign at the spike contributes maximally.
* **61D:** use spike amplitude 29/5 and physical tail height −7/4 at
  label 7, +7/4 at label 5, both of norm squared 38281/400.
  The third carrier uses spike amplitude 43/8 and tail height −3/5
  at label 6, of norm squared 141201/1600. This lower-amplitude
  carrier permits three labels at the same spike position. The
  extra tail dot with the old marked sparse banks is included,
  separately for the plus and minus original marked heads.

For donors with the same spike, use the distance-24 subcode; for
different spike sites, use the mother-code distance 16 because
their cubics are different. For old sparse banks, the complete
weight-16 support is orthogonal to the entire mother code. Thus the
previous marked and unmarked signs either force at least one body
mismatch, reducing the dot by at least two at their unit amplitude,
or are safe already without this saving (one of the marked cases
in 61D uses opposite tails). The same exact test accounts for the
ordinary and tail affine planes. All vectors are radially normalized
only after checking these norm-scaled inequalities.

Run python3 verify.py without -O. It pins and replays the sibling
transfer proof, constructs the DG generator, and checks every finite
algebraic angle, parity, distance and counting condition with exact
arithmetic. It requires the inherited Python/C++ verification
environment. An optional first positional argument selects a directory
holding the sibling sources. The calculation is a lower bound, not an
optimality or global-priority claim.
