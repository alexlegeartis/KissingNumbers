#!/usr/bin/env python3
"""One command to check, regenerate and re-verify everything.  Run this after any change.

    python update.py            # the fast path, about 20 minutes
    python update.py --full     # also the slow verifications, about two hours
    python update.py --check    # change nothing; just report whether it is all consistent
    python update.py --skip-paper --repository-only # maintain constructions without the paper

WHAT IT DOES, IN ORDER

  1. every package driver, so the shipped size tables and certificates are turned into claims
  2. audit.py, which recomputes each claim from those tables and cross-checks RESULTS.md
  3. audit.py --write-results, which regenerates RESULTS.md (including the status column)
  4. each package's make_readme_table.py, so the README tables come from the drivers
  4b. paper/mkvalues.py, so every bound the manuscript prints comes from RESULTS.md
  5. gpu/propagate.py, which installs a finished GPU class family into the package
     (skipped when ../gpu/ is absent, as it is in a clone of this repository)
  6. run_all.py, which re-runs every verification -- including the paper's checkers,
     factcheck.py and unsupported.py, which VERIFY its tables against RESULTS.md
     rather than regenerate them; the paper's tables have no generator

Steps 1-5 only move numbers that a driver computed; nothing here decides anything.  If step 2
reports a disagreement, that is the thing to look at -- RESULTS.md and the drivers have
diverged, and step 3 would paper over it.

HOW TO IMPROVE A BOUND.  RESULTS.md's status column says what would have to change:

    cosets     dimensions 36-37 and 73-87. Larger phase codes, nearest subsets,
               or compatible tail assignments remain open; no construction ceiling is claimed.
    heads      dimensions 25-31, 34-35, 38 and 40-46. Joint layer and family/head searches and partial families
               remain open. Complete-family populations retain their documented assumptions.
    codes      dimensions 32-33, 39, 61-68 and 95-96. Improve the explicit input codes,
               sparse support layers, compression maps or central sections.
    classes    dimensions 49-60 and 88-94. Better pairwise-disjoint class
               families improve the corresponding cap constructions. The October P48
               family is certified by dim49-63-p48-caps/verify_october.py.
    exhausted  dimension 18. These construction templates have
               their documented ceilings; other construction ideas remain possible.
    lp-exact   dimensions 69-71. The existing k-point programs are solved exactly;
               a different section or a larger k is needed to improve their bounds.

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
PAPER_DIR = os.path.join(HERE, 'paper')
V = os.path.join(HERE, 'verifications', 'improved')
FULL = '--full' in sys.argv
CHECK = '--check' in sys.argv
SKIP_PAPER = '--skip-paper' in sys.argv
AUDIT_ARGS = ['--repository-only'] if '--repository-only' in sys.argv else []

DRIVERS = [
    ('dims 36, 37', os.path.join(V, 'dim36-37-coxeter-todd'), 'final.py'),
    ('dims 41-46', os.path.join(V, 'dim41-46-p48-heads'), 'final.py'),
    ('dims 49-63', os.path.join(V, 'dim49-63-p48-caps'), 'final.py'),
    ('dims 73-80', os.path.join(V, 'dim73-95-gamma72-caps'), 'final73-80.py'),
    ('dims 81-95', os.path.join(V, 'dim73-95-gamma72-caps'), 'final81-95.py'),
    ('dims 46, 47', os.path.join(V, '..', 'recovered', 'dim46-47-p48-cross-sections'),
     'derive.py'),
]
DRIVERS += [(_name, os.path.join(V, _name), 'final.py') for _name in
            ('dim32-35-signed-layers', 'dim25-collective-deformation', 'dim31-dual-centres', 'dim38-joint-latitudes', 'dim39-deformed-orbits', 'dim40-affine-shells', 'dim29-30-midpoint-replacements', 'dim38-latitude-extension', 'dim39-root-replacements', 'dim26-27-midpoint-replacements', 'dim28-joint-leech-layers', 'dim61-68-95-96-code-layers', 'dim64-65-signed-parity', 'dim73-88-gamma72-cosets')]
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
r = run('audit.py', HERE, ['audit.py'] + AUDIT_ARGS)
if r is not None and r.returncode == 0:
    for l in r.stdout.split('\n'):
        if 'dims ' in l and any(s in l for s in ('classes', 'exhausted', 'external', 'lp-exact', 'heads')):
            print("     %s" % l.strip())

if not CHECK:
    print()
    print("3. regenerate RESULTS.md")
    run('audit.py --write-results', HERE, ['audit.py', '--write-results'] + AUDIT_ARGS)

    print()
    print("4. package README tables")
    for cwd, drv in TABLES:
        run(os.path.basename(cwd), cwd, [os.path.join('scripts', 'make_readme_table.py')])

    print()
    print("5. the paper's bounds: paper/values.tex from RESULTS.md")
    # Every claimed bound appears in kissing46.tex as \Kh{d}, defined once in values.tex.
    # The paper's PROSE and the addends of its count identities still have no generator, and
    # its checkers run in step 6 with the rest of run_all.py; only the bounds are regenerated
    # here.  Rebuild the PDF by hand afterwards if a bound moved -- see paper/README.md.
    if SKIP_PAPER:
        print('  skipped (--skip-paper): the manuscript remains its earlier snapshot')
    else:
        run('mkvalues.py', PAPER_DIR, ['mkvalues.py'], optional=True)

    print()
    print("6. install a finished GPU run, if one is waiting")
    # The workspace helper also rewrites manuscript tables, so skipping only
    # mkvalues.py would not keep the paper unchanged.
    if SKIP_PAPER:
        print('  skipped (--skip-paper): propagate.py also rewrites manuscript tables')
    else:
        run('propagate.py', os.path.join(ROOT, 'gpu'), ['propagate.py'], optional=True)
else:
    print()
    print("3-5 skipped (--check)")

print()
print("6. verifications%s" % (" (--full)" if FULL else ""))
run('run_all.py', HERE, ['run_all.py'] + (['--full'] if FULL else []) + AUDIT_ARGS +
    (['--skip-paper'] if SKIP_PAPER else []))

print()
print("=" * 92)
if fail:
    print("FAILED: %s" % ', '.join(fail))
else:
    print("EVERYTHING CONSISTENT AND VERIFIED")
print("=" * 92)
sys.exit(1 if fail else 0)
