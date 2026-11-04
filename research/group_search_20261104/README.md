# Group-representation discovery campaign, 4 November 2026

**Outcome: no new kissing-number record.** The 10D discovery control succeeds;
the 14D construction control is exactly reproduced. Higher-dimensional interfaces
and orbit-compressed checks work, but the tested Fourier families are far below
the current records. No software was installed, no changes were pushed, and no
Actions were enabled. All pre-existing repository files were preserved.

## What was actually implemented

`group_search.py` constructs even/odd subrepresentations of the finite Fourier and
quadratic-phase action on functions on F_p^r. It searches coordinate seeds, simple
character eigenspaces of sampled cyclic subgroups, and one-dimensional common
character eigenspaces of pairs of sampled group elements. It enumerates their
projective orbits, reconstructs **vector orbits with the actual phases**, and
selects compatible phase-rotated orbit unions by binary MILP.

This is a sampled representation/subgroup/character search, not a coordinate
insertion search. It is not an enumeration of all finite groups, all subgroups,
all characters, or all phase choices. The phase grid is multiples of pi/12,
limited to twelve offsets, and the admitted phase orders are 1,2,3,4,6,2p.
Orbit limits and floating-point rounding are explicit in the outputs. These can
miss seeds or configurations; negative results establish no global bound.

The distinction between projective lines and vector phases matters in practice:
canonicalizing each of the 40 lines independently and adding hexagons gives a
violation of about 0.57735. Reconstructing the true 80-vector orbit before adding
three phases yields the valid 240-vector component. A regression test checks it.

## Controls and results

| Experiment | Result | Scope |
|---|---:|---|
| 10D, even p=3,r=2; 120 sampled words, line cap 1000 | 510 | Discovered 40- and 45-line systems from generators; assembled 3x80+270 without input record coordinates |
| 14D, E7 orthogonal-frame construction | 1932 | Rebuilt Ganzhinov section 5.5 from its recipe, including an exact-cover search over 135 frames; **not** autonomous discovery of U(3,3) subgroup characters |
| 8D, odd p=3,r=2 | 240 | Numerical discovery/control |
| 10D, omit generator 0/1/4 | 96 / 510 / 174 | Three sampled subgroup-generator variants; no improvement |
| 28D, even p=3,r=3; 12 words, line cap 1600 | 6720 | Numerical 3x2240 orbit union; rank 28; all-pairs maximum cosine 0.5 within floating error |
| 26D, odd p=3,r=3; 12 words, line cap 1600 | no admitted orbit | Every sampled projective orbit exceeded the cap; **not** a nonexistence result |

The rank-one survey also used p=5,7,11,13,17 and both parity choices as listed in
`survey.log`, covering additional ambient dimensions 6,8,10,12,14,16. All results
are below records. There was no targeted search in real dimensions 5,7,9,11,13,15;
complex realification naturally prioritizes even dimensions in this first pass.
The separate `probe26.json` increases the vector cap to 20000 for one seed.

The solver's `upper_bound` is only for the numerically constructed finite block
graph; it is neither a rigorous mathematical upper bound nor a kissing-number
upper bound. A `best: 0` output means that this run admitted no block, not that
the dimension has zero kissing number.

## Exact independent verification

`verify.py` generates 10D using rational arithmetic in Q(sqrt(-3)); the seed of
the 270-vector orbit is the algebraic reconstruction of the numerically discovered
seed. It exports integer coordinates with diagonal quadratic form diag(I_9,3I_9).
Real coordinates are `C[:,j]*sqrt(weights[j]/norm_squared)`. Although this file
has 18 columns, its exact real rank is 10. Choosing any orthonormal basis of its
span realizes the code in R^10; no numerical rank estimate is needed.

The 14D certificate uses integer coordinates divided by sqrt(8). The verifier
checks distinctness, all norms, every pairwise inequality, and rank in exact
arithmetic. See `exact_verification.json`:

- 10D: 510 vectors, rank 10, 129795 pairs, norm 144, maximum weighted inner product 72.
- 14D: 1932 vectors, rank 14, 1865346 pairs, norm 8, maximum inner product 4.
- No violations. Certificate coordinate hashes are in the JSON.

`test_search.py` additionally checks rejection of duplicate points and a
norm-preserving integer mutation with a violated pair inequality, as well as
agreement between orbit-compressed and full cross-pair checks. Numerical search
NPZs are **not certificates**, even when their maximum cosine rounds to 0.5.

## Higher dimensions and scaling

The generator API accepts any unitary representation, including real orthogonal
ones and reducible representations. An external NPZ with a `generators` array of
shape (k,d,d) can be used with `--generators FILE --tag NAME`. The caller must
establish finiteness; the code checks unitarity and caps orbit enumeration.
No catalogue package or external solver was installed.

For transitive orbits A=Gv, B=Gw, unitary invariance gives
`max Re<a,b> = max_{b in B} Re<v,b>`.
A union of k scalar-phase copies needs only k representative rows against B,
rather than |A| rows. The implementation uses that identity for the compatibility
graph and still performs a separate full pair check on the selected code. It
never constructs the full group multiplication table or a dense |G| by |G| Gram
matrix. Double-coset enumeration, character tables and certified interval
arithmetic are **not yet implemented**.

The live Cohn table fetched on 2026-11-04 agrees with this repository's verified
values: **26D 201136; 28D 221988**. The 28D result is only 6720/221988 (~3.03%) of
its comparison. Baseline rows, URL and source hash are in `baselines.json`;
repository evidence is `common/source_review_20261029_27_31.json` and the
26D and 27–31D verification packages. The low-dimensional baselines remain
40,72,126,240,306,510,604,841,1154,1932,2564,4320 for dimensions 5–16.

### Next research decisions

1. Extend the representation catalogue beyond Fourier families, prioritizing
   low-degree sporadic-group and exceptional reflection-group representations
   in real dimensions 24–32, and reducible extensions in dimensions 25–31.
   Generate new seed vectors from subgroup-character projectors before selecting
   orbit unions; do not assume the Leech shell or a record equator is retained.
2. Replace random word pairs with subgroup classes and character projectors
   `(1/|H|) sum_h conjugate(nu(h)) rho(h)`, including higher-multiplicity spaces
   with continuous parameters. The present simple/common-eigenline sampling is
   only a first layer of this programme.
3. Compute cross-angle sets using stabilizers and double cosets, and solve
   continuous relative-phase constraints instead of relying on a fixed grid.
4. Permit splitting a dense full orbit into smaller subgroup orbits; rejecting
   a full orbit does not reject its useful subsets. The three subgroup trials
   above establish this interface, not an exhaustive strategy.
5. Require exact/interval certification and a renewed source comparison before
   any new-record claim. A larger numerical output alone is insufficient.

These are remaining research work, not completed features. There is no external
permission blocker. GAP/Sage and a representation catalogue were unavailable;
none was installed or requested. The mathematical objective of a new SOTA is
still open. In particular, the first higher-dimensional family is not presently
competitive, so merely increasing its sampling budget is not justified.

## Reproduce

Uses existing Python with NumPy, SciPy (MILP/HiGHS), SymPy; standard library otherwise.
From this directory:

```sh
OPENBLAS_NUM_THREADS=1 python verify.py
OPENBLAS_NUM_THREADS=1 python test_search.py
OPENBLAS_NUM_THREADS=1 python group_search.py --trials 120 --limit 1000
OPENBLAS_NUM_THREADS=1 python group_search.py --p 3 --rank 3 --parity 1 --trials 12 --limit 1600
OPENBLAS_NUM_THREADS=1 python group_search.py --p 3 --rank 3 --parity -1 --trials 12 --limit 1600
```

Search RNG seed is 20261104; exact seeds and generation formulas are in source.
Floating eigenspace ordering and MILP tie-breaking can vary across library versions;
exact certified cardinalities/inequalities do not depend on these numerical choices.

Sources: Ganzhinov, *Highly symmetric lines*, arXiv:2207.08266, sections 5.4–5.5,
https://arxiv.org/html/2207.08266 ; current baseline
https://cohn.mit.edu/kissing-numbers/ . The repository's `KNOWLEDGE.md` sections
63 and 127 were consulted to avoid repeating the closed 510 lifting/shell routes
and the 1932 ceiling in the norm-4 shell of E7+E7. Statements about those restricted
shells must not be generalized to all lattices or all configurations.

## Extension: complex code constructions

The follow-up branch in [CODE_CONSTRUCTIONS.md](CODE_CONSTRUCTIONS.md) now
reconstructs the DGLS ring lift from binary factors and discovers its coset
doubling. `code_search.py` and its exact receipts extend this campaign beyond
spherical group orbits; no new SOTA resulted.
