# Clone a dense owner onto multiple signed axes in 63D and 64D

| Dimension | GF(8) predecessor | New lower bound | Net gain |
|---:|---:|---:|---:|
| 63 | 184,427,081 | **184,427,691** | 610 |
| 64 | 335,796,250 | **335,797,303** | 1,053 |


This is an alternative to the preceding GF(8) pole construction, not an
addition to it. Keep the original cyclic dense code used in the verified
predecessors and all their sparse/signed layers. Select a list of dense
owners and a set of marked axes per owner. For each owner, remove it once,
and insert all its selected spikes. Retire the two original axial vectors
on each selected axis. The removed dense owners all existed in the original
configuration; in 63D each is in its published balanced octad section.

The vector for a word u on a marked axis a uses the sign word u xor e_a,
i.e. reverses the heavy sign only. Each ordinary magnitude is 2/sqrt(7);
the heavy magnitude is 14/sqrt(7). Thus squared norm is 64. Reversing the
heavy sign makes clones of the same dense owner on distinct axes compatible.
The verifier checks every new point pair with exact integer products in
units of 4/7. Against every other dense word, at least 16 parent sign
differences are guaranteed by the inherited cyclic code; the largest
possible product is 40 times 2/sqrt(7), strictly below 32.

Old four- and nine-support families are kept: both give the same bound
40 times 2/sqrt(7), and retained axes have bound 8 times 2/sqrt(7).
For a weight-sixteen signed vector, the only possible conflict is a
support containing the heavy coordinate with all 16 signs matching.
On the unrepaired classes the original sign words have even total parity,
as does every dense owner on each original support. Reversing the heavy
sign makes the new pattern odd on that support. Every exceptional
weight-sixteen class is checked directly using its pinned parity checks
and its affine shift.

For the 63D section, the heavy coordinate lies off the marked octad and
the owner has exactly four negative signs there. The new point therefore
lies in the same section. Reordering an unmodified standard class to
count that section does not change its even total sign parity.

Run `OPENBLAS_NUM_THREADS=1 python verify.py` with NumPy. The witness
files contain one row per removed dense owner, its unsigned binary sign
word and the bitmask of marked axes, both written as decimal integers.
