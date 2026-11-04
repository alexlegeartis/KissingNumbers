# Extending the campaign beyond simple group orbits

**Implemented outcome:** a second search branch reconstructs the published DGLS
(96,2^37,24) code from binary cyclic factors, Hensel lifting, a generalized Gray
map and discovered compatible translates. No new code or kissing record is claimed.
This is recovery within an explicitly supplied grammar of mathematical operations,
not autonomous invention of those operations.

## Source investigation and its limits

The requested 2001 work is *A Z8-linear lift of the binary Golay code and a nonlinear
binary (96,2^37,24)-code*, DOI 10.1109/18.923742. Its authors' university record
confirms the construction sequence and parameters:
https://experts.illinois.edu/en/publications/a-zsub8sub-linear-lift-of-the-binary-golay-code-and-a-nonlinear-b/ .
The original IEEE full text could not be fetched (HTTP 418 from the stamp endpoint;
the browser tool also could not open the document). The author's old publications
page returned 403. These were access failures, not approval denials; no bypass or
installation was attempted.

I read the primary follow-up by Duursma and Greferath, *Computing Symmetrized
Weight Enumerators for Lifted Quadratic Residue Codes*, arXiv:math/0302132,
https://arxiv.org/abs/math/0302132 . It explains why code descriptions and
weight-enumerator calculations can scale when enumerating all words cannot.
The repository's two independent exact parent verifiers provide a computational
check of the actual constructed code, separate from source metadata.

The follow-up also identifies an important next computational tool: partial,
character-weighted disjoint-support enumerators and transforms connecting ring
levels. A raw unsigned count of disjoint codeword pairs is **not** that enumerator;
coefficients can be negative. These recurrences have not been implemented here.

## What the new executable searches

`code_search.py` implements this sequence of typed operations:

1. Factor `1+x+...+x^(n-1)` over F2 and select half-degree factors.
2. Compute the unique coprime monic lift of each factor pair to Z4, Z8 and Z16,
   solving correction equations over F2 at each level. No published ring
   generator or lift coefficients are used as input.
3. Generate a free cyclic module and extend it by minus the coordinate sum.
4. Construct the affine Boolean Gray alphabet. Check **every symbol pair** for
   homogeneous-distance isometry before transferring any distance conclusion.
5. Search binary translation cosets of the local Gray alphabet; test repeated
   row masks using an additive lower bound on cross distances. Select compatible
   masks greedily. This selection is restricted and is not globally optimal.
6. Map binary words to signs. For length N and distance d, the resulting unit
   vectors have inner products at most `1-2d/N`; a kissing layer requires `4d>=N`.

The run explores n=7 and n=23, both half-degree factors, and ring moduli 4,8,16:
12 construction programs. Small modules of at most 2^18 words are exhaustively
checked using their homogeneous-weight histograms. Large unchecked programs are
explicitly marked `structural_candidate_only`.

## DGLS discovery control

After synthesis, one generated n=23,Z8 matrix exactly matches the repository's
published-code certificate. This comparison is **after** the factor/lift search.
The new generated input is passed to both existing parent verifiers:

- `verifications/improved/dim88-96-z8-golay/source/replay.py::verify_z8`
- `verifications/improved/dim88-96-z8-golay/source/verify_independent.py::parent_proof`

Their complete receipts agree: both check 2^24 Z4 words and 759 octad residual
cosets, each against 5489 small errors. They establish distance 24 for 2^36 Gray
images. The translation search finds two local cosets with certified cross
distance 24 over 24 rows, yielding **2^37 = 137438953472** distinct binary words.
No enumeration of 2^37 words is performed or claimed.

This control now goes earlier than the old repository replay: it reconstructs
the ring lift rather than taking the Z8 generator as a supplied input. It does
not replace the existing two independent verifiers or claim a new proof method.

## Other measured outputs

Both reciprocal n=7 factors give the following certified construction parameters:

| Binary length / ambient dimension | Parent cardinality | Parent distance | Kissing-oriented union size | Certified union distance |
|---|---:|---:|---:|---:|
| 16 | 256 | 6 | 256 | 6 |
| 32 | 4096 | 10 | 8192 | 8 |
| 64 | 65536 | 20 | 524288 | 16 |

These are generated controls, not novelty claims. In particular, 8192 in 32D is
far below the current 352800 baseline. The DGLS sign layer alone has 137438953472
points in 96D, below this repository's verified **142162362560** after additional
layers. This campaign found no improvement to either benchmark.

The distinction between parent and union distance is intentional. Preserving
distance 10 in the 32-bit parent admits one repeated-mask coset in this search;
the kissing objective only needs distance 8 and admits two. Construction search
must optimize the final geometric objective, not only isolated code parameters.

## General search architecture

The campaign now has two actual branches:

- **Representation branch:** unitary generators, subgroup-character seeds,
  orbit generation, phases, compatible orbit unions.
- **Ring-code branch:** binary factorization, chain-ring lifting, free modules,
  Gray embeddings, compatible translates and implicit sign layers.

Both output a construction recipe, cardinality, metric guarantee and verification
status. Compressed objects should remain compressed: a cardinality of 2^37 is a
reason to use algebraic certificates, not to allocate a giant point array.

Further operators should be added as separately tested modules:

- More finite rings and Galois rings; alternative residue/torsion codes and
  extension/gluing rules. Coprime Hensel lifting is unique for a fixed pair, so
  useful choices live in the residue factor, ring, extension and assembly—not
  repeated optimization of that same lift.
- Generalized concatenation with varying inner alphabets, code partitions and
  outer codes; nonconstant translate masks and switching.
- Symmetrized/partial weight-enumerator transforms, syndrome exclusions and
  automorphism reductions to certify distance without full enumeration.
- Geometric sections and additional signed-support layers, evaluated against
  the combined construction rather than only the dense sign layer.
- Group actions on code modules, coordinates, cosets and supports, alongside
  groups acting directly on spherical vectors.

These are **planned**, not implemented capabilities. The existing 88–96D package
already handles specific sections and extra layers of its fixed parent; this
campaign has not turned those into a general optimizer. Neither branch can yet
claim coverage of arbitrary constructions as sophisticated as DGLS. The useful
new guarantee is narrower: the implemented grammar can reconstruct this example
and transfer its exact metric proof through several construction stages.

## Reproduce and evidence

```sh
OPENBLAS_NUM_THREADS=1 python code_search.py
OPENBLAS_NUM_THREADS=1 python test_code_search.py
```

- `code_search_results.json`: all twelve recipes, statuses, small-code histograms,
  coset choices and the exact agreed DGLS parent receipt.
- `dgls_synthesized_generator.json`: generator and polynomial factors produced by
  synthesis, suitable as input to the independent existing verifier.
- `code_tests.log`: five passing tests, including an independent direct all-pairs
  binary distance check for the 256-word control, all Gray symbol pairs,
  incompatible-factor rejection, and objective-dependent translate selection.
- `code_search.log`: execution output. SymPy emits a harmless modular-factor
  ordering deprecation warning; it does not affect exact factorization checks.

No pre-existing repository file was modified. No push, Actions or software
installation. Remaining blockers are research scope and original-paper access,
not pending authorization.
