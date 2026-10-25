# Gamma72 nearest-coset subsets and repaired lifts

**Retained construction and history for dimensions 87–88.** The compact symmetry
family in [the cap package](../dim73-95-gamma72-caps/) now gives stronger bounds
there. This package still supplies the registered constructions in 73–86D.

Checkpoint: **25 October 2026**. The supplied verified constructions improve
dimensions **73–88**. The 88D construction independently passed all
11,640 pair families and has **6,231,205,108** points.

The driver `final.py` counts the strongest supplied verified construction in
each dimension. The current totals appear in [RESULTS.md](../../../RESULTS.md).
These are constructive bounds above the retained cap constructions. Cohn's
current table covers dimensions 1–48 and cannot establish priority here.

## Replay

```bash
python final.py
python verify73.py          # current 73D construction, a few seconds
python verify73.py --rebuild # also rebuild the short-vector enumerations
python verify.py
```

The replay requires NumPy, SymPy and threadpoolctl. It checks the complete
SHA-256 manifest, then copies the source into a temporary directory. It checks
the nearest-subset certificates, repaired lifts and all scaled cross-head
products. The full scaled replay can take about an hour, depending on the machine. All original reports and
coordinate files stay unchanged. `common/checkpoint.py` supplies the integrity
and temporary-copy helpers.

## Proof scope and inherited inputs

The parent is Nebe's extremal Gamma72 lattice, with minimum squared norm eight
and 6,218,175,600 minimal vectors. These parent facts are literature inputs.
The three seed sets have 6,008, 6,000 and 5,998 representatives of squared norm
22; division by two gives nearest-vector subsets of squared norm 11/2.

The nearestness proof uses 21 completed Leech-block interval enumerations.
The supplied lists and their completion metadata are retained. The replay
checks the rational LDL interval enclosures, list hashes and block populations,
and independently excludes every shorter triple needed by the construction.
It reconstructs the explicit seed heads and checks their norms, parity,
distinctness and angles. It does not rerun those 21 original enumerations;
`gamma72_leech_enum.cpp` is supplied for that additional step.

The initial lifts append compatible tails to nearest subsets. Later lifts
remove conflicts between exact automorphism images and tail assignments.
Integral isometry identities, distinct coset parities and exact tail ranks
are checked. Every retained cross-head inequality is independently regenerated;
search conflict caches are not used. Binary32 products are accepted only after
an integer absolute-term bound proves that every product and partial sum is
exactly representable. No construction ceiling is claimed.

## Files

`final.py` counts the verified totals. `verify.py` runs the isolated replays.
`UPDATE_PROVENANCE.json` records the supplied files and their source hashes.
`data/SHA256SUMS.txt` covers every installed file except itself.

`source/dim73-80/` contains the catalogue Gram and automorphism generators.
The remaining supplied files are under `source/dim61_96/`:

- `gamma72_coset_certificate.py` checks the seed nearestness and initial lifts.
- `gamma72_cosets.py` supplies rational LDL and coordinate helpers; its search
  entry points are not used by the replay.
- `gamma72_coset_certificate.json`, `gamma72_coset_heads.npz` contain the seed
  report and explicit representatives.
- `gamma72_leech_structure.npz`, `gamma72_leech_ldl.txt` give the block model
  and outward interval enclosures.
- The 21 `*.input.txt`, `*_r5_cert.bin` and matching JSON files give the block
  centres, complete short-vector lists and completion metadata.
- `gamma72_export_block*.bin`, `*.selected.bin` and
  `gamma72_head_*.triples.bin` reconstruct the explicit selected heads.
- `gamma72_leech_enum.cpp` is the interval enumerator; `gamma72_coset_screen.cpp`,
  `gamma72_coset_count.cpp`, `gamma72_coset_extract.cpp` preserve the original
  screening, counting and extraction code.
- `gamma72_repaired_verify.py`, `gamma72_repaired_lifts.npz` and its JSON check
  the repaired constructions in dimensions 75–79.
- `gamma72_scaled_verify_fast.py` checks the scaled constructions using exact
  binary32 products. `gamma72_scaled_verify.py` is the binary64 alternative.
- `gamma72_order91.npz`, `gamma72_scaled_seeds_*.npz`, `gamma72_scaled_*.npz`
  and their JSON reports give the exact group actions, tails, assignments
  and retained seed indices for dimensions 80–88.

The largest files are explicit proof data, including the completed block
enumerations and retained head lists. Finite search failures and unassigned
label families from the research workspace are not new kissing bounds.

The original 88D construction report is preserved as
`gamma72_scaled_16.pending.json`. The separate `gamma72_scaled_16.json`
records the completed independent replay. Its construction, seed, Gram and
verifier hashes match the installed files; `UPDATE_PROVENANCE.json` records
the replay and preserves the original report’s provenance.

## CPU conflict-cover continuation, 25 October 2026

The retained heads now use exact minimum vertex covers of the fixed search
conflict graphs. Maximum matching solves every bipartite component; CP-SAT
solves the components containing odd cycles. The 86D graph has 885,586 active
vertices, but only 777 lie in non-bipartite components. This reduces its cover
from 390,549 to 387,571 deletions and adds 11,912 points.

This pass improves 76D, 78D, 79D and every dimension from 80D through 86D.
The same search also improves the retained coset construction in 87D, where
the separate cap package remains stronger. All changed constructions passed
the original independent coordinate verifiers: 621,414,139 products for the
repaired lifts, 6,541,625,213 for the four-seed lifts and 26,222,187,021 for
the eight-seed lifts. The optimality statement concerns only these fixed
conflict graphs; other head sets and tail assignments remain open.

## Non-lattice improvement in dimension 73, 25 October 2026

The current bound is **6,218,191,606**, exceeding the completed lattice's
6,218,188,144 by **3,462**. The explicit lattice has determinant 5/2
and supplies 6,272 norm-22 representatives, whose two layers add 12,544
contacts to Gamma72. Its self-contained data and verifier are in
`source/dim73_surgery/lattice/`; `lattice/verify.py` checks the lattice,
and its README states all source facts and enumeration details.

The non-lattice construction retains those contacts. Its class contains
1,731 minimal lines u with |u.G.v| <= 2 within the class and
|u.G.x| <= 4 against all 6,272 representatives x. Remove each equatorial
pair +/-u and insert the four points

    (+/-sqrt(3)/2 * u, +/-sqrt(2)).

Each point has squared norm eight. Distinct class directions give products
at most 7/2; opposite tails above the same direction give four. Products
with retained equatorial vectors are at most 2 sqrt(3) < 4. Products with
the old coset contacts are at most sqrt(3)+sqrt(5) < 4. Thus each replaced
line gains two contacts, giving 6,218,188,144 + 2(1,731).

`source/dim73_surgery/verify_surgery.py` checks all class and mixed products
using int64 arithmetic. `surgery_class.npy` contains the class;
`surgery_search.json` records the CPU search and `surgery_verification.json`
its independent replay. `nonlattice_witness.npy` proves the configuration
cannot lie in a lattice: two retained equatorial vectors sum to a removed
u, while two caps sum to sqrt(3) u. These generate a dense subgroup on a line.
`verify73.py` checks the manifest and replays these checks in a temporary copy;
its `--rebuild` option also reruns `lattice/verify.py --rebuild`.

The new comparison is with the local certified lattice result, not a claim
that a comprehensive current literature table exists in dimension 73.
