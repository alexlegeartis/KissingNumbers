#!/usr/bin/env python3
"""Run every verification in this repository and report a single verdict.

    python run_all.py              # 93 scripts, including large exact certificate replays
    python run_all.py --full       # 112 scripts, budget about 10 hours
    python run_all.py --list       # just list what would run
    python run_all.py --changed    # do not re-run a job whose inputs have not moved
    python run_all.py --only dim31 # just the jobs whose label or command matches
    python run_all.py --skip-paper --repository-only # leave the manuscript snapshot alone

    python run_all.py --sota --quick # short replays for the results labelled SOTA in README.md
    python run_all.py --sota         # also replay the slow coset and binary-code checks

The SOTA profile omits historical constructions, searches and manuscript checks. Its
coverage report distinguishes a fresh replay, reuse of earlier checks, deferred work,
and inherited inputs. See SOTA_VERIFY.md for the remaining population and GPU-family gaps.

--changed.  A verification that re-runs unchanged work is mostly a waste of an hour, so each
job records what it read and the digest of every one of those files, in .run_all_cache.json,
and is not re-run while all of them and its own package tree are unchanged and the interpreter
and library versions are the same.  WHAT A JOB READS IS MEASURED, NOT DECLARED: these scripts
reach across package boundaries -- dim68-69 reads dim73-95's Gram, dim32-44-ers-audit reads
RESULTS.md, dim49-63 reads a coordinate file from outside the repository altogether -- so a
hand-written list would be wrong the first time anyone added a read.  Every child therefore runs
with an audit hook (sys.addaudithook, the "open" event) that logs each file it opens, and an
import is charged to the SOURCE its .pyc came from.  A job whose entry has no measured set --
because it passed before this existed -- is skipped only when NOTHING in the tree has changed.
Uncertainty means run.

The verdict says how many were not re-run and when they last passed; a run that skipped
anything never claims more than it checked.  --changed is opt-in, and it is not what a release
does: for that, and for anyone reproducing the repository, run the whole thing.

The minutes in the table below are BUDGETS, deliberately generous, so their sum overstates
the total -- the fast set is quoted from a measured run instead (2327 s over 59 scripts,
2026-09-16, with the machine to itself; 4495 s over the same 59 the afternoon before, and
1557 s over 55 on one laptop at the 2026-08-30 release).  RUN IT ALONE: the same 59
scripts took 10844 s on that afternoon while another job held a few gigabytes,
and dim29-30-frame-layer/verify.py 30 alone went from 174 s to 8518 s -- 49 times -- because
its blocked matmuls start swapping.  Nothing failed, but a wall-clock budget elsewhere in
the suite can, so a contended run is not a verdict on anything except the machine.

Two budgets were not generous but wrong.  scripts/verify_poles.py was budgeted at 0.5
minutes, which was right when the axis layers were the lattice-shell ones (248 points at
k = 23) and 19x low once they became rotated copies (93 074).  It is still the longest job
in the fast set, at 517 s measured against a budget of 18 minutes -- and it has been over
its own budget twice, at 9.4 and again at 10.5, each time because the layers under it grew.
dim31-sqrt3-layer/verify31.py was shipped at 3.0 and measured at 205 s alone and 282 s in
the suite, so it was over its budget on the day it was written; it went to 8.0, and on
2026-09-17 the package was superseded by dim31-frame-layer/verify.py 31, budgeted at 7.0
against 409 s measured alone.  A budget
that is only ever written once is a budget that describes the day it was written.

Nor is the TOTAL comparable across days: the machine moves under it, and so does the suite.
Compare the per-job times, not the sum.

Each script exits non-zero if any of its own checks fails; this collects the results.
Nothing here re-proves anything that the individual scripts do not -- it is a convenience,
and the place to look when something breaks after an edit.

Requirements: Python 3.9+, numpy and scipy for the verifications, and sympy for the
paper's formula check. The October coset replay also uses threadpoolctl. The optional older population replay additionally needs
python-flint and a C++17 g++ compiler. OR-Tools is imported only by CP-SAT models that are not
run from here -- the retired ones kept for the record, and
verifications/improved/dim73-95-gamma72-caps/scripts/k11_optimal.py, which settles dimension
83's partition as exactly optimal rather than merely best-found.  A job that cannot import one of
the packages named in OPTIONAL below is reported as skipped, and the skip says which package;
a job that cannot import anything else has FAILED -- above all a module of this repository's
own, which used to be reported as a skip inside a run that reported 0 FAILED.

MEMORY.  Two of the `--full` jobs are memory-hungry and are better run on their own:

    verifications/superseded/dim27-triple-partition/scripts/verify_exhaustive.py    (all
        20 108 045 530 pairs, blocked, but the blocks are large)
    verifications/improved/dim49-63-p48-caps/scripts/regenerate.py                (~4 GB)

On a machine with little free memory, `--full` can lose the whole run to the first of
those.  Run them individually if that happens; each is self-contained and exits non-zero on
failure like everything else here.
"""
import hashlib
import json
import os, re
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
V = os.path.join(HERE, 'verifications')


def _paper():
    """the paper's directory, inside this package or beside it -- whichever holds it

    It shipped INSIDE from 2026-08-30 so that its three checkers run in a clone instead of
    skipping; before that it was a sibling, and a checkout may still be arranged that way.
    Both are accepted, and a layout with neither leaves the three jobs to skip as they did."""
    for c in (os.path.join(HERE, 'paper'), os.path.join(HERE, os.pardir, 'paper')):
        if os.path.exists(os.path.join(c, 'kissing46.tex')):
            return c
    return os.path.join(HERE, 'paper')


PAPER = _paper()

# (label, directory, argv, minutes, in the fast set?)
JOBS = [
    ("audit: all 70 claims, mutual and external consistency",
     HERE, ['audit.py'], 0.5, True),
    ("common: extremal theta series and one-point distributions",
     os.path.join(HERE, 'common'), ['theta.py'], 0.1, True),
    ("common: the EXACT vertex, against three COUNTED Leech cross-sections",
     os.path.join(HERE, 'common'), ['exact_vertex.py'], 0.2, True),
    ("common: the k-point LP, validated on the Leech and P_48",
     os.path.join(HERE, 'common'), ['kpoint_lp.py'], 0.2, True),
    ("common: ADVERSARIAL -- every known truth inside its own LP bracket",
     os.path.join(HERE, 'common'), ['validate_lp.py'], 3.5, True),
    ("common: the cap construction's algebra and count, exactly",
     os.path.join(HERE, 'common'), ['capalgebra.py'], 0.1, True),

    # The paper ships in this repository, at paper/, so these three run for everybody.  They
    # were skipped in every clone until 2026-08-30, because the paper was kept outside on the
    # grounds that it cites the repository -- which is true and was not a reason to withhold
    # the 900-odd checks that tie the manuscript to the computations.  They still SKIP, rather
    # than fail, in a checkout that has no paper: verified by running this file from a copy of
    # the package with paper/ removed.
    ("paper: every figure in kissing46.tex, recomputed",
     PAPER, ['factcheck.py'], 0.1, True),
    ("paper: no figure in kissing46.tex is unsupported by the repository",
     PAPER, ['unsupported.py'], 0.1, True),
    ("paper: every formula in kissing46.tex, re-derived",
     PAPER, ['formulas.py'], 0.2, True),

    ("dim 25: the 197579-point lens-head configuration, exactly",
     os.path.join(V, 'improved', 'dim25-lens-heads'), ['verify.py'], 4.0, True),
    ("dim 18: the 8358-point bent-coset hexagon over the odd Barnes-Wall lattice, exactly",
     os.path.join(V, 'improved', 'dim18-bent-hexagon'), ['verify18.py'], 0.5, True),
    ("dim 25: the same, as an independent floating-point net over all 197579 points",
     os.path.join(V, 'improved', 'dim25-lens-heads'), ['fullcheck.py'], 2.0, True),
    ("dim 26: the 199806-point two-LAYER configuration, exactly",
     os.path.join(V, 'improved', 'dim26-27-iota-triangles'), ['verify26.py'], 1.0, True),
    ("dim 27: the 201557-point two-LAYER configuration, exactly",
     os.path.join(V, 'improved', 'dim26-27-iota-triangles'), ['verify27.py'], 2.0, True),
    ("dim 28: the 204896-point configuration with the norm-8 frame layer, exactly",
     os.path.join(V, 'improved', 'dim28-norm8-frame-layer'), ['verify28.py'], 1.0, True),
    ("dim 28: the direction weight and the axis are at their ceilings, exhaustively",
     os.path.join(V, 'improved', 'dim28-norm8-frame-layer'), ['ceilings.py'], 0.3, True),
    ("dim 29: the 209968-point configuration with the norm-8 frame layer, exactly",
     os.path.join(V, 'improved', 'dim29-30-frame-layer'), ['verify.py', '29'], 4.5, True),
    ("dim 30: the 221012-point configuration with the norm-8 frame layer, exactly",
     os.path.join(V, 'improved', 'dim29-30-frame-layer'), ['verify.py', '30'], 4.5, True),
    ("dims 29, 30: the direction weight and the axis are at their ceilings, exhaustively",
     os.path.join(V, 'improved', 'dim29-30-frame-layer'), ['ceilings.py', '29'], 0.3, True),
    ("dim 30 again: the same report at k = 6",
     os.path.join(V, 'improved', 'dim29-30-frame-layer'), ['ceilings.py', '30'], 0.3, True),
    ("dim 31: the 238662-point configuration with the norm-8 frame layer, exactly",
     os.path.join(V, 'improved', 'dim31-frame-layer'), ['verify.py', '31'], 7.0, True),
    ("dim 31: the direction weight and the axis are at their ceilings, exhaustively",
     os.path.join(V, 'improved', 'dim31-frame-layer'), ['ceilings.py', '31'], 0.5, True),
    ("superseded 27: the 200540-point triple-partition configuration, from its coordinate file",
     os.path.join(V, 'superseded', 'dim27-triple-partition'),
     ['scripts/verify_configuration.py'], 0.2, True),
    ("dim 38: the Leech cap construction at codimension 14",
     os.path.join(V, 'improved', 'dim38-leech-large-codimension'),
     ['scripts/verify.py'], 0.2, True),
    ("dim 38: the axis rotation RE-DERIVED, descent then Cayley (a search: slow)",
     os.path.join(V, 'improved', 'dim38-leech-large-codimension'),
     ['scripts/axis_rotation.py'], 4.0, False),
    ("dims 38, 39: the Edel-Rains-Sloane chain and its codes",
     os.path.join(V, 'improved', 'dim39-ers-constant-weight'),
     ['scripts/verify.py'], 0.2, True),
    ("coxeter-todd: dims 36, 37, exact phase orbits and two latitude layers",
     os.path.join(V, 'improved', 'dim36-37-coxeter-todd'), ['verify.py'], 3.0, True),
    ("p48-heads: dims 41-46, partial 43D and coordinated 45D, inherited populations",
     os.path.join(V, 'improved', 'dim41-46-p48-heads'), ['verify.py'], 3.0, True),
    ("p48-heads: older 41/42/44 population certificate, fresh replay in a temporary copy",
     os.path.join(V, 'improved', 'dim41-46-p48-heads'),
     ['verify.py', '--replay-populations'], 20.0, False),
    ("p48-heads: recovered 43D seed, fresh populations and exact angles (2878155 points)",
     os.path.join(V, 'improved', 'dim41-46-p48-heads'),
     ['certify43.py'], 20.0, False),
    ("dims 46, 47: cross-sections of P_48, recovered not claimed",
     os.path.join(V, 'recovered', 'dim46-47-p48-cross-sections'), ['derive.py'], 0.1, True),
    ("october-p48: historical and refilled classes, all isometries and disjointness",
     os.path.join(V, 'improved', 'dim49-63-p48-caps'), ['verify_october.py'], 1.0, True),
    ("Z8 Golay: dimensions 88-96, two exact compressed-construction proofs",
     os.path.join(V, 'improved', 'dim88-96-z8-golay'), ['verify.py'], 5.0, True),
    ("october-codes: compression, support layers and central sections",
     os.path.join(V, 'improved', 'dim61-68-95-96-code-layers'), ['verify.py'], 5.0, True),
    ("dims 26-31: residual completions, full independent exact replay",
     os.path.join(V, 'improved', 'dim26-31-residual-completions'), ['verify.py'], 30.0, True),
    ("dims 27-31: joint refills, complete exact geometry and direct integer checks",
     os.path.join(V, 'improved', 'dim27-31-joint-refills'), ['verify.py'], 3.0, True),
    ("dims 26, 28, 29, 30: asymmetric exchanges, complete exact replay",
     os.path.join(V, 'improved', 'dim26-30-asymmetric-exchanges'), ['verify.py'], 10.0, True),
    ("dim27: quadratic normal module, complete independent all-pairs replay",
     os.path.join(V, 'improved', 'dim27-quadratic-normal-module'), ['verify.py'], 5.0, True),
    ("dim31: normal-module exchange, complete independent all-pairs replay",
     os.path.join(V, 'improved', 'dim31-normal-module-exchange'), ['verify.py'], 5.0, True),
    ("dims 26, 27: shared midpoint replacements, exact independent replay",
     os.path.join(V, 'improved', 'dim26-27-midpoint-replacements'), ['verify.py'], 1.0, True),
    ('dims 32-35: signed layers and coordinated heads, two independent exact proofs',
     os.path.join(V, 'improved', 'dim32-35-signed-layers'), ['verify.py'], 1.0, True),
    ('dim25-collective-deformation: independent exact replay', os.path.join(V, 'improved', 'dim25-collective-deformation'), ['verify.py'], 1.0, True),
    ('dim31-dual-centres: independent exact replay', os.path.join(V, 'improved', 'dim31-dual-centres'), ['verify.py'], 5.0, True),
    ('dim38-joint-latitudes: independent exact replay', os.path.join(V, 'improved', 'dim38-joint-latitudes'), ['verify.py'], 1.0, True),
    ('dim39-deformed-orbits: independent exact replay', os.path.join(V, 'improved', 'dim39-deformed-orbits'), ['verify.py'], 1.0, True),
    ('dim40-affine-shells: independent exact replay', os.path.join(V, 'improved', 'dim40-affine-shells'), ['verify.py'], 5.0, True),
    ("dims 29, 30: midpoint replacements, exact independent replay",
     os.path.join(V, 'improved', 'dim29-30-midpoint-replacements'), ['verify.py'], 2.0, True),
    ("dim 38: latitude extension and exact rational axis rotation",
     os.path.join(V, 'improved', 'dim38-latitude-extension'), ['verify.py'], 0.2, True),
    ("dim 39: neutral support changes and exact root replacements",
     os.path.join(V, 'improved', 'dim39-root-replacements'), ['verify.py'], 0.2, True),
    ("dims 73, 77-79, 82-85, 87: norm-ten addition and exact existence proofs",
     os.path.join(V, 'improved', 'dim73-88-gamma72-cosets'), ['verify_norm10.py'], 0.3, True),
    ("dim 28: joint Leech layers, exact independent replay",
     os.path.join(V, 'improved', 'dim28-joint-leech-layers'), ['verify.py'], 1.0, True),
    ("dims 64, 65: signed support augmentation and BCH owner caps",
     os.path.join(V, 'improved', 'dim64-65-signed-parity'), ['verify.py'], 1.0, True),
    ("october-codes-full: all 2^36 words of the 73-coordinate binary code",
     os.path.join(V, 'improved', 'dim61-68-95-96-code-layers'), ['verify.py', '--full-codes'], 7.0, False),
    ("dim 73: completed lattice plus non-lattice cap surgery",
     os.path.join(V, 'improved', 'dim73-88-gamma72-cosets'), ['verify73.py'], 0.1, True),
    ("october-cosets: Gamma72 nearest subsets and all repaired/scaled cross products",
     os.path.join(V, 'improved', 'dim73-88-gamma72-cosets'), ['verify.py'], 60.0, True),
    ("dims 49-63: the P_48 lattice, class and arithmetic",
     os.path.join(V, 'improved', 'dim49-63-p48-caps'), ['verify.py'], 0.1, True),
    ("dims 49-63: the result table",
     os.path.join(V, 'improved', 'dim49-63-p48-caps'), ['final.py'], 0.1, True),
    ("dims 49-63: regenerate and re-verify all 1400 classes",
     os.path.join(V, 'improved', 'dim49-63-p48-caps'), ['scripts/regenerate.py'], 8.0, False),
    ("dims 62, 63: the Edel-Rains-Sloane sign code (level 0)",
     os.path.join(V, 'improved', 'dim62-63-ers-chain'), ['scripts/verify.py'], 0.1, True),
    ("dims 62, 63: the FULL chain, levels 1 and 2 built, tau(8) and tau(16) reproduced",
     os.path.join(V, 'improved', 'dim62-63-ers-chain'), ['scripts/verify_chain.py'], 0.6, True),
    ("dims 57-63: the cap-direction line systems, exactly (lambda(k) = tau(k)/2)",
     os.path.join(V, 'improved', 'dim49-63-p48-caps'), ['scripts/verify_capdirs.py'], 0.3, True),
    ("dims 49-63: the cap construction's absolute ceiling, against level 0",
     os.path.join(V, 'improved', 'dim49-63-p48-caps'), ['scripts/ceiling.py'], 0.1, True),
    ("dims 70, 71: cross-sections of Gamma_72",
     os.path.join(V, 'improved', 'dim70-71-gamma72-cross-sections'), ['derive.py'], 0.2, True),
    ("dims 73-95: CALIBRATION against Cohn's table, dims 25-31",
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'), ['scripts/calibrate.py'], 2.0, True),
    ("dims 73-95: zero-sum triple partitions of A_1 ... E_8",
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'), ['scripts/triples.py'], 0.3, True),
    ("k-point LP: TIGHT on the Leech through k = 4, against counted embeddings",
     os.path.join(V, 'improved', 'dim68-69-gamma72-cross-sections'),
     ['scripts/leech_truth.py'], 0.1, True),
    ("dim 69: the A_3 triple at the Hermite floor, exactly",
     os.path.join(V, 'improved', 'dim68-69-gamma72-cross-sections'),
     ['scripts/verify69.py'], 0.1, True),
    ("dim 68: the programme's EXACT optimum, 361275480, by the dual vertex",
     os.path.join(V, 'improved', 'dim68-69-gamma72-cross-sections'),
     ['scripts/verify68.py'], 5.5, True),
    ("dims 68, 69, 70: every LP certificate re-derived by a second arithmetic path",
     os.path.join(V, 'improved', 'dim68-69-gamma72-cross-sections'),
     ['scripts/recheck_cert.py'], 1.0, True),
    ("dim 69: the bound does not depend on which basis presents the Gram",
     os.path.join(V, 'improved', 'dim68-69-gamma72-cross-sections'),
     ['scripts/presentation_invariance.py'], 1.5, True),
    ("dim 68: the four-point LP certificate, both bases",
     os.path.join(V, 'improved', 'dim68-69-gamma72-cross-sections'),
     ['scripts/k4d4.py'], 15.0, False),
    ("dims 81-95: cap-direction configurations and partitions, exactly",
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'),
     ['scripts/verify_parts.py'], 2.0, True),
    ("dims 84-85: the Kappa cap sets K_12 and K_13 RE-DERIVED from the Golay code",
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'),
     ['scripts/kappa_caps.py'], 4.0, True),
    ("dims 84-85: their axis layers RE-DERIVED, one plane rotation per block (a search)",
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'),
     ['scripts/kappa_poles.py'], 20.0, False),
    ("dims 84, 88, 90: the perfect triple partitions RE-DERIVED from an order-3 isometry",
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'),
     ['scripts/omega_parts.py', '12', '16', '18'], 6.0, True),
    ("dims 92, 94: the same, at rank 20 and 22 (a search; the k = 22 Gram is 1.2e9 pairs)",
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'),
     ['scripts/omega_parts.py', '20', '22'], 60.0, False),
    ("dims 89, 91, 93: the perfect triple partitions of the ODD k RE-DERIVED (a search)",
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'),
     ['scripts/oddparts.py', '17', '19', '21'], 8.0, True),
    ("dim 95: the same at k = 23, whose triple table is 115 million rows (a search)",
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'),
     ['scripts/oddparts.py', '23'], 60.0, False),
    ("the identity of every shipped cap set, as a digest",
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'),
     ['scripts/capsha.py'], 0.2, True),
    ("dims 81-95: all 32000 classes pairwise disjoint (needs scripts/classes)",
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'),
     ['scripts/verify_all32000.py'], 2.0, False),
    ("dim 95: the 950 extension classes, full Gram each (needs scripts/classes)",
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'),
     ['scripts/verify_ext.py'], 1.5, False),
    ("dims 74-80: the root-lattice axis layers, every legal shell, exactly",
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'),
     ['scripts/poles18.py'], 0.2, True),
    ("dims 75-95: the axis rotations RE-DERIVED, Cayley over each metric (a search)",
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'),
     ['scripts/axis_rotate.py'], 25.0, False),
    ("dims 75-95: the axis layers, exactly",
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'),
     ['scripts/verify_poles.py'], 18.0, True),
    ("dims 91-93: the axis layers from the RECORD configuration, RE-DERIVED (a search)",
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'),
     ['scripts/axis_record.py'], 120.0, False),
    ("dim 83: the axis layer over Z[sqrt2], RE-DERIVED (a search)",
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'),
     ['scripts/axis_sqrt2.py'], 30.0, False),
    ("dims 85-90: the axis layers from the PUBLISHED record configurations, RE-DERIVED",
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'),
     ['scripts/axis_cohn.py'], 90.0, False),
    ("dims 85-93: the record-configuration axis layers, exactly",
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'),
     ['scripts/verify_record.py'], 12.0, True),
    ("dim 83: the k = 11 configuration, exactly in Z[sqrt2]",
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'),
     ['scripts/verify11.py'], 0.5, True),
    ("dims 73-95: the 120 Gamma_72 classes, exactly",
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'),
     ['scripts/verify_classes.py', 'data/disjoint_classes.npz'], 0.5, True),
    ("dims 81-95: the 32000-class family REBUILDS and matches its shipped sizes",
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'),
     ['scripts/verify_rebuild.py', '60'], 1.5, True),
    ("dims 86-94: the GPU class family's certificate, without a GPU",
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'),
     ['scripts/gpu/check_cert.py'], 0.1, True),
    ("dims 87-94: compact symmetry classes, cap directions, poles and counts, on CPU",
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'),
     ['scripts/verify_symmetry.py', '--complete'], 2.0, True),
    ("dims 73-80: the result table",
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'), ['final73-80.py'], 0.5, True),
    ("dims 81-95: the result table",
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'), ['final81-95.py'], 0.1, True),

    ("superseded 44, 45: the LP brackets, and Sun-Wang inside them",
     os.path.join(V, 'superseded', 'dim44-45-p48-cross-sections'), ['brackets.py'], 40.0, False),

    ("closed 9-19: every record is maximal (needs Cohn's data set)",
     os.path.join(V, 'closed', 'dim09-19-record-maximality'), ['sweep.py'], 1.5, True),
    ("closed 22, 23: Lambda_21/22/23 are maximal",
     os.path.join(V, 'closed', 'dim22-23-maximal-cross-sections'), ['verify.py'], 0.1, True),
    ("closed 17-24: the layer identity, as a structure inside the Leech lattice",
     os.path.join(V, 'closed', 'dim17-24-layer-identity'), ['layers.py'], 0.6, True),
    ("closed 17-24: the flat layer exactly, and dims 20 and 21 rebuilt from the Golay code",
     os.path.join(V, 'closed', 'dim17-24-layer-identity'), ['flats.py'], 1.8, True),
    ("closed 17-24: the deep-hole radius rule, and A(16,8,6) = 16",
     os.path.join(V, 'closed', 'dim17-24-layer-identity'), ['equator.py'], 0.4, True),
    ("closed: the class problem, 425 bound exactly, by two routes",
     os.path.join(V, 'closed', 'class-problem-upper-bound'), ['certificate.py'], 0.1, True),
    ("closed 25: the ceiling of the lens-head construction, 196560 + 971 + alpha + 3",
     os.path.join(V, 'closed', 'dim25-construction-ceiling'), ['verify.py'], 17.0, False),
    ("closed 32-44: the Edel-Rains-Sloane audit",
     os.path.join(V, 'closed', 'dim32-44-ers-audit'), ['PIPELINE.py'], 0.1, True),
    ("closed 96: the cap construction at k = 24, and its Leech triples exactly",
     os.path.join(V, 'closed', 'dim96-ers-takeover'), ['cap96.py'], 0.1, True),
    ("closed 96: Edel-Rains-Sloane takes over, 12886999232 against 6480558568",
     os.path.join(V, 'closed', 'dim96-ers-takeover'), ['ers96.py'], 0.1, True),
    ("closed: the ERS level-0 floor in EVERY claimed dimension (62, 63 use it)",
     os.path.join(V, 'closed', 'dim96-ers-takeover'), ['ers_sweep.py'], 0.1, True),
    ("closed: the FULL ERS chain in every claimed dimension, and each claim's exposure",
     os.path.join(V, 'closed', 'dim96-ers-takeover'), ['ers_exposure.py'], 0.2, True),
    ("closed: why the table stops at 96 -- what 96 and a direct sum already imply",
     os.path.join(V, 'closed', 'dim96-ers-takeover'), ['above96.py'], 0.2, True),
    ("closed 19: the Cayley graph structure",
     os.path.join(V, 'closed', 'dim17-23-cohn-li-mechanism'),
     ['scripts/verify19.py'], 0.1, False),
    ("closed 17: alpha = 192 exactly, by a six-number Delsarte dual",
     os.path.join(V, 'closed', 'dim17-layered-family'),
     ['scripts/theta17.py'], 0.4, True),
    ("closed 17: the frame measured from the record (needs Cohn's data set)",
     os.path.join(V, 'closed', 'dim17-layered-family'),
     ['scripts/frame17.py'], 1.5, False),
]


# Jobs needed for the current README records, without rediscovering constructions or
# replaying their superseded versions. Do not infer geometric verification from a count
# driver, a manifest, or the GPU run's receipt.
def _sota_jobs(quick=False):
    wanted = {
        ('common', ('capalgebra.py',)),
        ('dim27-31-joint-refills', ('verify.py',)),
        ('dim26-30-asymmetric-exchanges', ('verify.py',)),
        ('dim27-quadratic-normal-module', ('verify.py',)),
        ('dim31-normal-module-exchange', ('verify.py',)),
        ('dim18-bent-hexagon', ('verify18.py',)),
        ('dim32-35-signed-layers', ('verify.py',)),
        ('dim25-collective-deformation', ('verify.py',)),
        ('dim31-dual-centres', ('verify.py',)),
        ('dim40-affine-shells', ('verify.py',)),

        ('dim26-31-residual-completions', ('verify.py',)),
        ('dim26-27-midpoint-replacements', ('verify.py',)),
        ('dim28-joint-leech-layers', ('verify.py',)),
        ('dim29-30-midpoint-replacements', ('verify.py',)),
        ('dim38-joint-latitudes', ('verify.py',)),
        ('dim39-deformed-orbits', ('verify.py',)),
        ('dim73-88-gamma72-cosets', ('verify_norm10.py',)),
        ('dim36-37-coxeter-todd', ('verify.py',)),
        ('dim41-46-p48-heads', ('verify.py',)),
        ('dim49-63-p48-caps', ('verify_october.py',)),
        ('dim49-63-p48-caps', ('final.py',)),
        ('dim49-63-p48-caps', ('scripts/verify_capdirs.py',)),
        ('dim61-68-95-96-code-layers', ('verify.py',)),
        ('dim88-96-z8-golay', ('verify.py',)),
        ('dim64-65-signed-parity', ('verify.py',)),
        ('dim68-69-gamma72-cross-sections', ('scripts/verify69.py',)),
        ('dim68-69-gamma72-cross-sections', ('scripts/recheck_cert.py',)),
        ('dim70-71-gamma72-cross-sections', ('derive.py',)),
        ('dim73-88-gamma72-cosets', ('verify73.py',)),
        ('dim73-95-gamma72-caps', ('scripts/verify_symmetry.py', '--complete')),
        ('dim73-95-gamma72-caps', ('final73-80.py',)),
        ('dim73-95-gamma72-caps', ('final81-95.py',)),
    }
    if not quick:
        wanted.remove(('dim61-68-95-96-code-layers', ('verify.py',)))
        wanted.add(('dim61-68-95-96-code-layers', ('verify.py', '--full-codes')))
        wanted.add(('dim73-88-gamma72-cosets', ('verify.py',)))
    current_packages = {package for _, package in _sota_records().values()}
    wanted = {job for job in wanted if job[0] == 'common' or job[0] in current_packages}
    jobs = []
    found = set()
    for lab, directory, argv, minutes, fast in JOBS:
        selector = (os.path.basename(directory), tuple(argv))
        if selector in wanted:
            found.add(selector)
            jobs.append((lab, directory, argv, minutes, fast))
    if found != wanted:
        raise ValueError('SOTA profile has no job for %r' % sorted(wanted - found))
    return [(JOBS[0][0], HERE, ['audit.py', '--repository-only'], 0.5, True)] + jobs


def _sota_records():
    """The current README labels are the scope; validate their values and packages.

    This deliberately does not assert that a literature search proves global priority.
    Adding a new labelled record without a verification route fails closed.
    """
    with open(os.path.join(HERE, 'RESULTS.md'), encoding='utf-8') as fh:
        registry = {}
        for line in fh:
            cells = line.split('|')
            if len(cells) < 7 or not cells[1].strip().isdigit():
                continue
            dimension = int(cells[1])
            value = int(''.join(c for c in cells[3] if c.isdigit()))
            package = re.search(r'\(verifications/improved/([^/)]+)/?\)', line)
            if not package or dimension in registry:
                raise ValueError('ambiguous registry row: %s' % line.strip())
            registry[dimension] = (value, package.group(1))
    records = {}
    with open(os.path.join(HERE, 'README.md'), encoding='utf-8') as fh:
        for line in fh:
            cells = line.split('|')
            if len(cells) < 11 or not cells[1].strip().isdigit() or cells[5].strip() != '**SOTA**':
                continue
            dimension = int(cells[1])
            value = int(''.join(c for c in cells[4] if c.isdigit()))
            if dimension in records or registry[dimension][0] != value:
                raise ValueError('SOTA/registry mismatch in dimension %d' % dimension)
            records[dimension] = registry[dimension]
    if not records:
        raise ValueError('README.md contains no recognised SOTA rows')
    return records


def _sota_coverage(jobs, results=None):
    """Report scope per dimension, including steps absent from the selected profile."""
    groups = [
        ('dim27-31-joint-refills', [27, 28, 29, 30, 31], [('verify.py',)], 'all norms and all unordered point pairs; bounded exact integer BLAS and separate direct-int64 checks of every inserted point'),
        ('dim26-30-asymmetric-exchanges', [26, 28, 29, 30], [('verify.py',)], 'complete explicit integer codes; exact opposite-pair and unpaired-point checks cover every unordered pair'),
        ('dim27-quadratic-normal-module', [27], [('verify.py',)], 'complete quadratic coordinates; all 5491467600 distinct antipodal-line pairs checked exactly without baseline assumptions'),
        ('dim31-normal-module-exchange', [31], [('verify.py',)], 'complete explicit integer coordinates; all 7881966235 distinct antipodal-line pairs checked without baseline assumptions'),
        ('dim18-bent-hexagon', [18], [('verify18.py',)], ''),
        ('dim32-35-signed-layers', [32, 33, 34, 35], [('verify.py',)], 'two independent exact representations; full dense-code enumeration and every exceptional pair type'),
        ('dim25-collective-deformation', [25], [('verify.py',)], ''),
        ('dim31-dual-centres', [31], [('verify.py',)], ''),
        ('dim40-affine-shells', [40], [('verify.py',)], ''),

        ('dim26-31-residual-completions', [26,27,28,29,30,31], [('verify.py',)], 'full independent coordinate replays and self-contained shell lemma; exact quadratic fields in 30D/31D'),
        ('dim26-27-midpoint-replacements', [26, 27], [('verify.py',)], ''),
        ('dim28-joint-leech-layers', [28], [('verify.py',)], ''),
        ('dim29-30-midpoint-replacements', [29, 30], [('verify.py',)], ''),
        ('dim36-37-coxeter-todd', [36, 37], [('verify.py',)], ''),
        ('dim38-joint-latitudes', [38], [('verify.py',)], ''),
        ('dim39-deformed-orbits', [39], [('verify.py',)], ''),
        ('dim41-46-p48-heads', list(range(41, 47)), [('verify.py',)],
         'inherited body populations; the latest large shells are not all freshly enumerated'),
        ('dim49-63-p48-caps', list(range(49, 61)),
         [('verify_october.py',), ('scripts/verify_capdirs.py',)], ''),
        ('dim61-68-95-96-code-layers', [61, 62, 63, 95, 96], [('verify.py',)], ''),
        ('dim64-65-signed-parity', [64, 65], [('verify.py',)], ''),
        ('dim61-68-95-96-code-layers', [66, 67, 68], [('verify.py',)],
         'the default replay reads the saved [73,36,16] distance report; --sota also enumerates all words'),
        ('dim68-69-gamma72-cross-sections', [69],
         [('scripts/verify69.py',), ('scripts/recheck_cert.py',)],
         'geometry replay plus an independent LP-certificate replay; recheck_cert also covers dimensions 68 and 70'),
        ('dim70-71-gamma72-cross-sections', [70, 71], [('derive.py',)], ''),
        ('dim73-88-gamma72-cosets', [73], [('verify73.py',)],
         'completed half-coset minimum, cap surgery, all refilled norm-ten heads in two exact arithmetics and non-lattice witness; published Gamma72 parent facts'),
        ('dim73-88-gamma72-cosets', [77, 78, 79, 82, 83, 84, 85, 87], [('verify_norm10.py',)],
         'exact finite-graph existence proofs; published Gamma72 minimum and all-shell design theorem; no full selected head list'),
        ('dim73-88-gamma72-cosets', [74, 75, 76, 80, 81, 86], [('verify.py',)],
         'seed nearestness inherits 21 completed block-enumeration lists and their completion metadata'),
        ('dim88-96-z8-golay', list(range(88, 97)), [('verify.py',)],
         'published nonlinear parent; complete finite distance proof, both support representations and exact section counts; no full 2^37-word enumeration'),
        ('dim73-95-gamma72-caps', list(range(88, 95)),
         [('scripts/verify_symmetry.py', '--complete')],
         'all retained class geometry and disjointness, directions, partitions, poles and counts; CPU replay'),
    ]
    records = _sota_records()
    # Retained historical routes may share dimensions with newer packages.
    # Coverage follows the current README package, and excludes surpassed rows.
    groups = [(package, [d for d in dimensions
                         if d in records and records[d][1] == package], commands, scope)
              for package, dimensions, commands, scope in groups]
    groups = [group for group in groups if group[1]]
    mapped = {}
    for package, dimensions, commands, scope in groups:
        for dimension in dimensions:
            if dimension in mapped:
                raise ValueError('two verification routes for dimension %d' % dimension)
            mapped[dimension] = package
    if set(records) != set(mapped):
        raise ValueError('SOTA coverage differs from README: missing %r, extra %r' % (
            sorted(set(records) - set(mapped)), sorted(set(mapped) - set(records))))
    if any(package != mapped[d] for d, (_, package) in records.items()):
        raise ValueError('SOTA package differs from its verification route')
    selected = {(os.path.basename(d), tuple(argv)): lab for lab, d, argv, _, _ in jobs}
    if results is not None:
        verdicts = {lab: verdict for verdict, lab, _ in results}
        audit_label = JOBS[0][0]
        common_label = JOBS[5][0]
    print('\nVerification scope for %d results labelled SOTA in README.md:' % len(records))
    for package, dimensions, commands, scope in groups:
        states = []
        for command in commands:
            lab = selected.get((package, command))
            if lab is None and command == ('verify.py',):
                lab = selected.get((package, ('verify.py', '--full-codes')))
            states.append('DEFERRED' if lab is None else
                          ('SELECTED' if results is None else verdicts.get(lab, 'DEFERRED')))
        if results is not None:
            states += [verdicts.get(audit_label, 'DEFERRED')]
            if package == 'dim49-63-p48-caps':
                states += [verdicts.get(common_label, 'DEFERRED')]
        state = next((s for s in ('FAIL', 'SKIP', 'DEFERRED', 'SAME') if s in states),
                     'SELECTED' if results is None else 'PASS')
        ds = ','.join(map(str, dimensions))
        print('  %-24s %-8s %s' % (ds, state, package))
        if scope:
            print('      %s' % scope)
    print('PASS covers the replay described above, with published parent lattice/code facts as inputs.')
    print('SAME reuses earlier checks after input hashing. DEFERRED and GAP are not fresh verification.')
    print('This profile does not check whether a later publication has surpassed a README record.')


# Third-party packages a run is allowed to be without.  A job that cannot import one of
# these has not failed, it has not run; anything else that fails to import HAS failed.
# This used to read `'ModuleNotFoundError' in out`, which also swallowed a module this
# repository owns: rename ers_exposure.py, leave an importer behind, and the job reported
# SKIP inside a run that reported 0 FAILED.  Name the module and check it against a list.
#
# sympy is used by paper/formulas.py in the fast set. python-flint is needed only
# by the optional original population replay in the full set. The repository also
# uses ortools, cvxpy, torch and triton,
# but only in scripts nothing here runs: the retired CP-SAT models and the GPU builders.
# Adding a name that no job can raise would read as coverage without being it, so a name
# belongs here when a job that needs it is registered, and not before.
OPTIONAL = ('sympy', 'flint')


def _absent_optional(out):
    """The optional package this output says is missing, or None."""
    for m in re.finditer(r"No module named '([^']+)'", out):
        top = m.group(1).split('.')[0]
        if top in OPTIONAL:
            return top
    return None



# ============================================================================================
# --changed: do not re-run a job whose inputs have not moved since it last passed
# ============================================================================================
# WHAT A JOB DEPENDS ON IS MEASURED, NOT DECLARED.  These scripts read across package
# boundaries -- dim68-69 reads dim73-95's Gram matrix, dim32-44-ers-audit reads RESULTS.md,
# dim49-63 reads a coordinate file from OUTSIDE the repository -- so a hand-written dependency
# list would be wrong the first time someone added a read.  Instead every child is run with
# PYTHONPATH pointing at a sitecustomize that installs sys.addaudithook and logs the path of
# every file the job opens.  That measured set, plus a digest of the job's own package tree
# (which catches files ADDED to it, not only files changed), is the key.
#
# A job with no measured set yet -- because it last passed before this mode existed -- may be
# skipped only when the digest of the WHOLE tree is unchanged.  Uncertainty always means run.

CACHE = os.path.join(HERE, '.run_all_cache.json')
_HOOK = r"""
import os, sys
_log = os.environ.get('RUNALL_OPENLOG')
if _log:
    try:
        _f = open(_log, 'a', buffering=1, encoding='utf-8', errors='replace')

        def _runall_hook(event, args):
            if event == 'open' and args and isinstance(args[0], str):
                try:
                    _f.write(os.path.abspath(args[0]) + '\n')
                except Exception:
                    pass
        sys.addaudithook(_runall_hook)
    except Exception:
        pass
"""
_SKIPDIR = ('__pycache__', '.git')


def _sha(path):
    h = hashlib.sha256()
    try:
        with open(path, 'rb') as fh:
            for blk in iter(lambda: fh.read(1 << 20), b''):
                h.update(blk)
    except OSError:
        return None
    return h.hexdigest()


def _tree_digest(root):
    """one digest over every file under root, names included, caches excluded"""
    h = hashlib.sha256()
    for dp, dn, fn in os.walk(root):
        dn[:] = sorted(d for d in dn if d not in _SKIPDIR)
        for f in sorted(fn):
            if f.endswith('.pyc') or f == os.path.basename(CACHE):
                continue
            p = os.path.join(dp, f)
            h.update(os.path.relpath(p, root).replace(os.sep, '/').encode('utf-8'))
            h.update((_sha(p) or 'gone').encode('ascii'))
    return h.hexdigest()


def _env_fingerprint():
    out = ['py%d.%d.%d' % sys.version_info[:3]]
    for m in ('numpy', 'scipy', 'sympy', 'ortools'):
        try:
            out.append('%s%s' % (m, __import__(m).__version__))
        except Exception:                                          # noqa: BLE001
            out.append('%s-' % m)
    return ' '.join(out)


def _load_cache():
    try:
        with open(CACHE, encoding='utf-8') as fh:
            c = json.load(fh)
        return c if isinstance(c, dict) else {}
    except Exception:                                              # noqa: BLE001
        return {}


def _save_cache(c):
    try:
        with open(CACHE, 'w', encoding='utf-8') as fh:
            json.dump(c, fh, indent=1, sort_keys=True)
    except OSError:
        pass


def _key(d, argv):
    return '%s | %s' % (os.path.relpath(d, HERE).replace(os.sep, '/'), ' '.join(argv))


def _unchanged(ent, d, envfp, gdig, selfcache):
    """(True, note) if this job's inputs have not moved since it passed"""
    if not ent or ent.get('status') != 'PASS':
        return False, ''
    if ent.get('env') != envfp:
        return False, ''
    deps = ent.get('deps')
    if deps is None:
        return (ent.get('tree') == gdig), 'tree'
    if d not in selfcache:
        selfcache[d] = _tree_digest(d)
    if ent.get('self') != selfcache[d]:
        return False, ''
    for rel, dig in deps.items():
        p = rel if os.path.isabs(rel) else os.path.join(HERE, rel)
        if _sha(p) != dig:
            return False, ''
    return True, '%d inputs' % len(deps)


def _hookdir():
    d = tempfile.mkdtemp(prefix='runall-hook-')
    with open(os.path.join(d, 'sitecustomize.py'), 'w', encoding='utf-8') as fh:
        fh.write(_HOOK)
    return d


def _traced_env(hookdir, logfile):
    e = dict(os.environ)
    e['RUNALL_OPENLOG'] = logfile
    pp = e.get('PYTHONPATH')
    e['PYTHONPATH'] = hookdir + (os.pathsep + pp if pp else '')
    return e


_PREFIXES = tuple(os.path.normcase(os.path.abspath(p))
                  for p in (sys.prefix, sys.base_prefix) if p)


def _source_of(p):
    """map <dir>/__pycache__/<mod>.cpython-39.pyc back to <dir>/<mod>.py, else p"""
    d, f = os.path.split(p)
    if os.path.basename(d) != '__pycache__' or not f.endswith('.pyc'):
        return p
    src = os.path.join(os.path.dirname(d), f.split('.')[0] + '.py')
    return src if os.path.isfile(src) else None


def _read_opens(logfile, limit=6000):
    """{path relative to HERE (or absolute if outside): digest}, or None if too many

    An import is recorded as the .pyc it actually read; that is mapped back to the source, so
    editing a module a job imports from ANOTHER package invalidates the entry.  Files under the
    interpreter and site-packages are dropped -- the environment fingerprint covers those, and
    they would swamp the list.
    """
    try:
        with open(logfile, encoding='utf-8', errors='replace') as fh:
            raw = {ln.strip() for ln in fh if ln.strip()}
    except OSError:
        return None
    out = {}
    for p in raw:
        p = _source_of(p)
        if not p or not os.path.isfile(p):
            continue
        ap = os.path.abspath(p)
        nc = os.path.normcase(ap)
        if nc.startswith(_PREFIXES) or 'site-packages' in nc or '__pycache__' in nc:
            continue
        if ap == CACHE or ap == os.path.abspath(logfile):
            continue
        try:
            rel = os.path.relpath(ap, HERE)
        except ValueError:
            rel = ap
        rel = ap if rel.startswith('..') else rel.replace(os.sep, '/')
        dg = _sha(ap)
        if dg:
            out[rel] = dg
        if len(out) > limit:
            return None
    return out


def main():
    full = '--full' in sys.argv
    sota = '--sota' in sys.argv
    quick = '--quick' in sys.argv
    if quick and not sota:
        print('--quick requires --sota; the ordinary suite is unchanged')
        return 2
    changed_only = '--changed' in sys.argv
    only = None
    if '--only' in sys.argv:
        i = sys.argv.index('--only')
        if i + 1 >= len(sys.argv):
            print("--only needs a substring to match against the label or the command")
            return 2
        only = sys.argv[i + 1].lower()
    jobs = _sota_jobs(quick) if sota else [j for j in JOBS if j[4] or full]
    if '--skip-paper' in sys.argv:
        jobs = [j for j in jobs if not j[0].startswith('paper:')]
    if '--repository-only' in sys.argv:
        jobs = [(lab, directory, argv + ['--repository-only'] if argv == ['audit.py'] else argv,
                 minutes, fast) for lab, directory, argv, minutes, fast in jobs]
    if only:
        jobs = [j for j in jobs if only in j[0].lower() or only in ' '.join(j[2]).lower()]
        if not jobs:
            print("--only %r matches no job; --list shows them all" % only)
            return 2
    if '--list' in sys.argv:
        for lab, d, argv, mins, fast in jobs:
            print("  %-58s ~%4.1f min   %s" % (lab, mins, ' '.join(argv)))
        print("")
        print("  %d scripts, about %.0f minutes" % (len(jobs), sum(j[3] for j in jobs)))
        if sota:
            _sota_coverage(jobs)
        return 0

    if sota:
        _sota_coverage(jobs)

    print("running %d scripts, roughly %.0f minutes%s%s"
          % (len(jobs), sum(j[3] for j in jobs), "" if full else "  (--full adds more)",
             "  [--changed: unchanged jobs are not re-run]" if changed_only else ""))
    print("=" * 92)
    results = []
    skipnote = {}
    cache = _load_cache()
    envfp = _env_fingerprint()
    selfcache = {}
    gdig = _tree_digest(HERE) if changed_only else None
    hookdir = _hookdir()
    logfile = os.path.join(hookdir, 'opens.log')
    same = []
    for lab, d, argv, mins, fast in jobs:
        t0 = time.time()
        sys.stdout.write("  %-62s " % lab[:62])
        sys.stdout.flush()
        if not os.path.isdir(d):
            print("%-5s %5s" % ("SKIP", "-"))
            results.append(("SKIP", lab, "directory not present: %s" % d))
            continue
        k = _key(d, argv)
        if changed_only:
            ok, note = _unchanged(cache.get(k), d, envfp, gdig, selfcache)
            if ok:
                print("%-5s %5s" % ("SAME", "-"))
                same.append((lab, cache[k].get('when', '?'), note))
                results.append(("SAME", lab, ""))
                continue
        try:
            if os.path.exists(logfile):
                os.remove(logfile)
        except OSError:
            pass
        try:
            p = subprocess.run([sys.executable] + argv, cwd=d, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, env=_traced_env(hookdir, logfile))
            out = p.stdout.decode('utf-8', 'replace')
            rc = p.returncode
        except Exception as e:                                    # noqa: BLE001
            out, rc = str(e), 99
        el = time.time() - t0
        if rc == 0 and 'SKIP --' in out:
            # the script ran, found its input absent and checked nothing
            verdict = "SKIP"
        elif rc == 0:
            verdict = "PASS"
        elif _absent_optional(out):
            verdict = "SKIP"
            skipnote[lab] = "needs %s" % _absent_optional(out)
        else:
            verdict = "FAIL"
        print("%-5s %5.0fs" % (verdict, el))
        results.append((verdict, lab, out))
        if verdict == 'PASS':
            if d not in selfcache:
                selfcache[d] = _tree_digest(d)
            cache[k] = {'status': 'PASS', 'when': time.strftime('%Y-%m-%d %H:%M'),
                        'env': envfp, 'self': selfcache[d], 'deps': _read_opens(logfile),
                        'tree': None, 'seconds': round(el, 1)}
        else:
            cache.pop(k, None)
        # after EVERY job, not at the end: these runs get interrupted, and a run that banks
        # nothing when it is stopped at job 60 of 67 is the thing --changed exists to avoid
        _save_cache(cache)

    print("=" * 92)
    _save_cache(cache)
    bad = [r for r in results if r[0] == 'FAIL']
    skipped = [r for r in results if r[0] == 'SKIP']
    for v, lab, out in bad:
        print("")
        print("FAILED: %s" % lab)
        print("\n".join(out.rstrip().split("\n")[-25:]))
    for v, lab, out in skipped:
        print("skipped: %s%s"
              % (lab, "  -- %s" % skipnote[lab] if lab in skipnote else ""))
    if same:
        print("")
        print("NOT RE-RUN -- inputs unchanged since the run that passed them (--changed):")
        for lab, when, note in same:
            print("   %-62s  passed %s%s"
                  % (lab[:62], when, "  (%s)" % note if note else ""))
    print("")
    print("%d passed, %d unchanged, %d skipped, %d FAILED"
          % (len(results) - len(bad) - len(skipped) - len(same), len(same),
             len(skipped), len(bad)))
    if same and not bad:
        print("   -- that verdict covers %d of %d jobs; the other %d passed earlier and "
              "nothing they read has moved" % (len(results) - len(same), len(results), len(same)))
    if sota:
        _sota_coverage(jobs, results)
    return 1 if bad else 0


if __name__ == '__main__':
    raise SystemExit(main())

