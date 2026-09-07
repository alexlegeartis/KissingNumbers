#!/usr/bin/env python3
"""Run every verification in this repository and report a single verdict.

    python run_all.py              # the fast set: 54 scripts, 26 minutes measured
    python run_all.py --full       # also the slow ones: 71 scripts, budget about 9 hours
    python run_all.py --list       # just list what would run

The minutes in the table below are BUDGETS, deliberately generous, so their sum overstates
the total -- the fast set is quoted from a measured run instead (1557 s on one laptop,
2026-08-30, the release).  One of them was not generous but wrong: scripts/verify_poles.py was
budgeted at 0.5 minutes, which was right when the axis layers were the lattice-shell ones
(248 points at k = 23) and 19x low once they became rotated copies (93 074).  It is still the longest job
in the fast set, at 517 s measured against a budget of 18 minutes -- and it has been over its
own budget twice, at 9.4 and again at 10.5, each time because the layers under it grew.  A
budget that is only ever written once is a budget that describes the day it was written.

Nor is the TOTAL comparable across days: the machine moves under it, and so does the suite.
Compare the per-job times, not the sum.

Each script exits non-zero if any of its own checks fails; this collects the results.
Nothing here re-proves anything that the individual scripts do not -- it is a convenience,
and the place to look when something breaks after an edit.

Requirements: Python 3.8+, numpy and scipy for the verifications, and sympy for the
paper's formula check.  Nothing else: OR-Tools is imported only by CP-SAT models that are not
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
import os, re
import subprocess
import sys
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
    ("audit: all 48 claims, mutual and external consistency",
     HERE, ['audit.py'], 0.1, True),
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

    ("dim 25: the 197569-point lens-head configuration, exactly",
     os.path.join(V, 'improved', 'dim25-lens-heads'), ['verify.py'], 4.0, True),
    ("dim 25: the same, as an independent floating-point net over all 197569 points",
     os.path.join(V, 'improved', 'dim25-lens-heads'), ['fullcheck.py'], 2.0, True),
    ("dim 26: the 199632-point two-triangle configuration, exactly",
     os.path.join(V, 'improved', 'dim26-27-iota-triangles'), ['verify26.py'], 1.0, True),
    ("dim 27: the 201010-point two-triangle configuration, exactly",
     os.path.join(V, 'improved', 'dim26-27-iota-triangles'), ['verify27.py'], 2.0, True),
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
    ("dims 46, 47: cross-sections of P_48, recovered not claimed",
     os.path.join(V, 'recovered', 'dim46-47-p48-cross-sections'), ['derive.py'], 0.1, True),
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
     os.path.join(V, 'improved', 'dim73-95-gamma72-caps'), ['scripts/calibrate.py'], 1.0, True),
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
    ("closed: the class problem, 425 bound exactly, by two routes",
     os.path.join(V, 'closed', 'class-problem-upper-bound'), ['certificate.py'], 0.1, True),
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
     ['scripts/verify19.py'], 1.0, False),
    ("closed 17: alpha = 192 exactly, by a six-number Delsarte dual",
     os.path.join(V, 'closed', 'dim17-layered-family'),
     ['scripts/theta17.py'], 0.4, True),
    ("closed 17: the frame measured from the record (needs Cohn's data set)",
     os.path.join(V, 'closed', 'dim17-layered-family'),
     ['scripts/frame17.py'], 1.5, False),
]


# Third-party packages a run is allowed to be without.  A job that cannot import one of
# these has not failed, it has not run; anything else that fails to import HAS failed.
# This used to read `'ModuleNotFoundError' in out`, which also swallowed a module this
# repository owns: rename ers_exposure.py, leave an importer behind, and the job reported
# SKIP inside a run that reported 0 FAILED.  Name the module and check it against a list.
#
# sympy is the whole list, because it is the only one any registered job imports -- paper/
# formulas.py, in the fast set.  The repository does use ortools, cvxpy, torch and triton,
# but only in scripts nothing here runs: the retired CP-SAT models and the GPU builders.
# Adding a name that no job can raise would read as coverage without being it, so a name
# belongs here when a job that needs it is registered, and not before.
OPTIONAL = ('sympy',)


def _absent_optional(out):
    """The optional package this output says is missing, or None."""
    for m in re.finditer(r"No module named '([^']+)'", out):
        top = m.group(1).split('.')[0]
        if top in OPTIONAL:
            return top
    return None


def main():
    full = '--full' in sys.argv
    jobs = [j for j in JOBS if j[4] or full]
    if '--list' in sys.argv:
        for lab, d, argv, mins, fast in jobs:
            print("  %-58s ~%4.1f min   %s" % (lab, mins, ' '.join(argv)))
        print("")
        print("  %d scripts, about %.0f minutes" % (len(jobs), sum(j[3] for j in jobs)))
        return 0

    print("running %d scripts, roughly %.0f minutes%s"
          % (len(jobs), sum(j[3] for j in jobs), "" if full else "  (--full adds more)"))
    print("=" * 92)
    results = []
    skipnote = {}
    for lab, d, argv, mins, fast in jobs:
        t0 = time.time()
        sys.stdout.write("  %-62s " % lab[:62])
        sys.stdout.flush()
        if not os.path.isdir(d):
            print("%-5s %5s" % ("SKIP", "-"))
            results.append(("SKIP", lab, "directory not present: %s" % d))
            continue
        try:
            p = subprocess.run([sys.executable] + argv, cwd=d, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
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

    print("=" * 92)
    bad = [r for r in results if r[0] == 'FAIL']
    skipped = [r for r in results if r[0] == 'SKIP']
    for v, lab, out in bad:
        print("")
        print("FAILED: %s" % lab)
        print("\n".join(out.rstrip().split("\n")[-25:]))
    for v, lab, out in skipped:
        print("skipped: %s%s"
              % (lab, "  -- %s" % skipnote[lab] if lab in skipnote else ""))
    print("")
    print("%d passed, %d skipped, %d FAILED"
          % (len(results) - len(bad) - len(skipped), len(skipped), len(bad)))
    return 1 if bad else 0


if __name__ == '__main__':
    raise SystemExit(main())
