#!/usr/bin/env python3
"""One command to check, regenerate and re-verify everything.  Run this after any change.

    python update.py            # the fast path, about 20 minutes
    python update.py --full     # also the slow verifications, about two hours
    python update.py --check    # change nothing; just report whether it is all consistent

WHAT IT DOES, IN ORDER

  1. every package driver, so the shipped size tables and certificates are turned into claims
  2. audit.py, which recomputes each claim from those tables and cross-checks RESULTS.md
  3. audit.py --write-results, which regenerates RESULTS.md (including the status column)
  4. each package's make_readme_table.py, so the README tables come from the drivers
  5. gpu/propagate.py, which installs a finished GPU class family into the package
     (skipped when ../gpu/ is absent, as it is in a clone of this repository)
  6. run_all.py, which re-runs every verification -- including the paper's checkers,
     factcheck.py and unsupported.py, which VERIFY its tables against RESULTS.md
     rather than regenerate them; the paper's tables have no generator

Steps 1-5 only move numbers that a driver computed; nothing here decides anything.  If step 2
reports a disagreement, that is the thing to look at -- RESULTS.md and the drivers have
diverged, and step 3 would paper over it.

HOW TO IMPROVE A BOUND.  RESULTS.md's status column says what would have to change:

    classes    36 of the 47, dimensions 49-61 and 73-95.  Find a better family of
               pairwise-disjoint classes and they all improve together.  P_48:
               dim49-63-p48-caps/scripts/regenerate.py, whose family is automorphism images of
               one base class -- more classes, or images chosen to overlap less.  Gamma_72:
               the GPU builder, which certifies a family by a seed.  Its scripts are
               shipped at verifications/improved/dim73-95-gamma72-caps/scripts/gpu/; the
               working tree that RUNS them, and the integrate.py that installs a
               finished run, sit in ../gpu/ beside this package and are not
               part of it.  Before spending compute, read
               dim49-63-p48-caps/scripts/ceiling.py: the construction has an ABSOLUTE ceiling
               per dimension, and in 62 and 63 that ceiling is already below what a single
               sign code gives, which is why those two are not in this group.
    external   dimensions 39, 62, 63 and 96.  Each is a function of a PUBLISHED TABLE and needs
               no new mathematics: recheck Brouwer's constant-weight tables for 39, and
               codetables.de for [62,27,16] or [63,28,16], either of which would double its
               dimension.  Today both stand at d = 15, one short.  Dimension 96 rests on
               A(96,24) >= 2^33 from [96,33,24], two thirds of its total, and on the grid
               maps [24,9,12]_4 and RS[6,4,3]_16; only the last is built here.
    exhausted  dimensions 25, 26, 27, 38.  The idea is finished; only a different construction
               (in 25: a plateau traverse of the lens-head search moves H by one sphere at a time)
               will move these.
    lp-exact   dimensions 68-71.  NOT "a Gram the search has not reached" -- the Gram is
               forced.  A rank-k sublattice of Gamma_72 has det >= (8/gamma_k)^k by Hermite,
               the bound falls as the determinant rises, and k = 1..4 all sit at the floor
               (A_2, A_3, D_4), each exhibited.  What is left is (a) the BASIS, which changes
               the numeric dual and so the certificate by up to 0.7% -- scan_k3.py and
               scan_k4.py sweep the attaining Grams -- and (b) a larger k, which needs the
               construction rewritten in array form and extrapolates below the floor at k = 5.

Nothing here knows about tau(n) itself: no matching upper bound is known in any of these
dimensions, so no claim is known to be optimal.  What IS known is when a claim cannot be
improved by its own construction, and audit.py checks two such ceilings: the level-0
Edel-Rains-Sloane floor under every claim (step 4b), and, where two constructions reach the
same dimension, that the better one is the one claimed (step 4a).
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
V = os.path.join(HERE, 'verifications', 'improved')
FULL = '--full' in sys.argv
CHECK = '--check' in sys.argv

DRIVERS = [
    ('dims 49-63', os.path.join(V, 'dim49-63-p48-caps'), 'final.py'),
    ('dims 73-80', os.path.join(V, 'dim73-95-gamma72-caps'), 'final73-80.py'),
    ('dims 81-95', os.path.join(V, 'dim73-95-gamma72-caps'), 'final81-95.py'),
    ('dims 46, 47', os.path.join(V, '..', 'recovered', 'dim46-47-p48-cross-sections'),
     'derive.py'),
]
TABLES = [(os.path.join(V, 'dim49-63-p48-caps'), 'final.py'),
          (os.path.join(V, 'dim73-95-gamma72-caps'), 'final81-95.py')]

fail = []


def run(label, cwd, args, optional=False):
    t0 = time.time()
    if not os.path.isdir(cwd):
        if not optional:
            print("  %-46s MISSING %s" % (label, cwd))
            fail.append(label)
        return None
    r = subprocess.run([sys.executable] + args, cwd=cwd, capture_output=True, text=True)
    ok = r.returncode == 0
    print("  %-46s %s  %5.0fs" % (label, "ok  " if ok else "FAIL", time.time() - t0))
    if not ok:
        fail.append(label)
        tail = [l for l in (r.stdout + r.stderr).strip().split('\n') if l.strip()][-6:]
        for l in tail:
            print("        %s" % l[:120])
    return r


print(__doc__.split('WHAT IT DOES')[0])
print("1. package drivers")
for label, cwd, script in DRIVERS:
    run(label, cwd, [script], optional=True)

print()
print("2. audit: recompute every claim and cross-check RESULTS.md")
r = run('audit.py', HERE, ['audit.py'])
if r is not None and r.returncode == 0:
    for l in r.stdout.split('\n'):
        if 'dims ' in l and any(s in l for s in ('classes', 'exhausted', 'external', 'lp-exact')):
            print("     %s" % l.strip())

if not CHECK:
    print()
    print("3. regenerate RESULTS.md")
    run('audit.py --write-results', HERE, ['audit.py', '--write-results'])

    print()
    print("4. package README tables")
    for cwd, drv in TABLES:
        run(os.path.basename(cwd), cwd, [os.path.join('scripts', 'make_readme_table.py')])

    print()
    print("5. install a finished GPU run, if one is waiting")
    # NOT 'the paper': the paper's tables have no generator, and its checkers run in
    # step 6 with the rest of run_all.py.  This step installs a GPU class family.
    run('propagate.py', os.path.join(ROOT, 'gpu'), ['propagate.py'], optional=True)
else:
    print()
    print("3-5 skipped (--check)")

print()
print("6. verifications%s" % (" (--full)" if FULL else ""))
run('run_all.py', HERE, ['run_all.py'] + (['--full'] if FULL else []))

print()
print("=" * 92)
if fail:
    print("FAILED: %s" % ', '.join(fail))
else:
    print("EVERYTHING CONSISTENT AND VERIFIED")
print("=" * 92)
sys.exit(1 if fail else 0)
