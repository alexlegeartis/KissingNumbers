# Independent review of the compatible sign families

`verify_independent.py` is a second implementation of the family checks,
63D section count, 61D/62D compression inequalities, and 65D cap transfer. It
imports no discovery scripts or other verification modules. Its construction
inputs are the local `data/global64.json` and `data/short61_62.json`, plus the
pinned parent code archive and dual-octad list named in `data/source_inputs.json`.

Run from the repository root, respecting the available CPU allocation:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verifications/improved/dim61-67-compatible-sign-families/verify_independent.py
```

The public `check()` function returns a JSON-compatible report. The CLI prints
that report and writes no files. Neither entry point rewrites the recorded
receipt, `independent_verification.json`. Report fields are deterministic;
elapsed time and machine-specific paths are omitted. The entry point rejects
Python's `-O` mode, and the actual checks use explicit exceptions.

## Replayed results

| Dimension | Independently recounted points |
|---:|---:|
|61|88,975,146|
|62|122,539,020|
|63|184,393,988|
|64|335,711,104|
|65|335,779,330|

The global family has 1,940 new supports. The reviewer checks 61,687,150
new/old and new/new support pairs, finding 15,520 old/new and 303 new/new
ten-overlap repairs. The shortened family has 937 new supports. Its 9,502,117
changed support pairs include 4,046 old/new and 99 new/new repairs.

Every repair is independently reconstructed from the intersection of the two
explicit dual codes: a common weight-eight check must lie inside the support
intersection and have opposite affine phases. Every sign dual is checked to be
RM(1,4), and its kernel is reconstructed by eliminating columns over the binary
field. The all-one parity check and distinct syndrome columns forbid nonzero
kernel words of weight one, two, or three; the resulting full sign class has
2,048 words and minimum distance four.

The section reviewer checks all 72 supplied dual octads, verifies the rank-seven
projection onto each, and chooses a best section without reading the production
section result. It then explicitly enumerates 7,643,136 affine sign words for
the chosen octad `1284508757262368`. Its terms are:

| Layer | Points |
|---|---:|
|Dense|146,800,640|
|Retained old weight-sixteen supports|35,503,776|
|New weight-sixteen supports|1,986,560|
|Four-support layer|102,900|
|Poles|112|
|**Total**|**184,393,988**|

The compression check uses rational upper and lower enclosures for square
roots, independently of the main verifier's quadratic-field comparison. It
also checks the dense projection ranks and resulting cardinalities, all 736
61D dense pair cases and 36 hexagon pair types, the dense/four-support
interactions, the affine-plane quadruple design, and the 65D transfer bounds.

## Choices requiring an explicit statement

For 63D, a previously untouched old Hamming sign code may be reordered so that
its intersection with the chosen octad is a dual hyperplane. This keeps all
2,048 original signs before taking the section. These supports meet every new
support in at most eight positions, so their reordered signs preserve all
cross-support bounds. The section count uses that freedom.

For 65D, the dense layer is the separately verified BCH radical used by the
owner construction. It is not the 64D dense code. Every interaction with the
changed sparse layers follows from amplitude bounds, and none of the 68,224
owners belongs to the weight-sixteen layer. The transfer therefore preserves
the owner construction.

## Scope and other exact replays

This reviewer does not duplicate the exhaustive parent-code enumeration, the
65D BCH owner replay, or the full 66D/67D baseline checks. Those are performed
by `check_family.cpp`, `replay_parents.py`, `check_expanded.cpp`, and
`baseline67.py` through the main package verifier. The retained parent support
pairs in the expanded dimensions are checked directly, so no minimum-distance
claim about every word of the 73-bit mother code is needed for those layers.

The 67D proper-subcode addition has its own independent implementation,
`verify_subcode67_independent.py`, which directly multiplies the affine sign
pairs for every exceptional support overlap. Its separate receipt records the
49-class, 4,624-point gain. No mathematical gap was found in either reviewed
construction. This is an exact verification of the supplied witnesses and
does not assert their optimality.
