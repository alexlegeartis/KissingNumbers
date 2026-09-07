#!/usr/bin/env python3
"""Run the whole verification suite.

    python scripts/run_all.py            the two checks that need no numpy
    python scripts/run_all.py --full     also the exhaustive all-pairs check
    python scripts/run_all.py --quick    only the checker that proves the theorem

Order:
  1. verify_configuration.py   proves the theorem from the data file alone
  2. verify_paper_claims.py    checks the paper's finite claims and the bridge
  3. selftest.py               negative controls, with --full
  4. verify_exhaustive.py      all 20108045530 pairs, with --full (needs numpy)
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))


def run(script, label):
    print()
    print("#" * 78)
    print(f"# {label}")
    print("#" * 78, flush=True)
    t = time.time()
    code = subprocess.call([sys.executable, os.path.join(HERE, script)])
    print(f"\n({script} finished in {time.time() - t:.0f}s, exit {code})", flush=True)
    return code


def main(argv):
    quick = "--quick" in argv
    full = "--full" in argv
    steps = [("verify_configuration.py", "1. the theorem, from the data file alone")]
    if not quick:
        steps.append(("verify_paper_claims.py", "2. the paper's finite claims, and the bridge"))
    if full:
        steps.append(("selftest.py", "3. negative controls"))
        steps.append(("verify_exhaustive.py", "4. all pairs, assuming nothing (needs numpy)"))

    t0 = time.time()
    failures = [label for script, label in steps if run(script, label) != 0]
    print()
    print("=" * 78)
    if failures:
        print(f"{len(failures)} STEP(S) FAILED:")
        for f in failures:
            print(f"  - {f}")
    else:
        print(f"All {len(steps)} steps passed in {time.time() - t0:.0f}s.")
        print()
        print("data/dimension27_200540.txt is a valid kissing configuration of 200540")
        print("points in R^27, so K(27) >= 200540.")
    print("=" * 78)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
