#!/usr/bin/env python3
"""Why the table stops at 96, measured rather than asserted.

    python above96.py

Dimension 96 is the last row of every table in this repository, and the reason usually
given -- that the layered construction stops winning there -- is the reason 96 belongs to
Edel-Rains-Sloane, not the reason 97 is absent.  The actual reason is that K is
superadditive: an n-dimensional and an m-dimensional 60-degree code, placed on orthogonal
subspaces, meet at inner product 0, so

    K(n + m) >= K(n) + K(m) ,

and every value above 96 is therefore already implied by dimension 96 together with the
published low-dimensional records.  A row for dimension 100 reading
K(96) + K(4) = 12 886 999 232 + 24 would be arithmetic on this paper's own dimension 96,
not a result -- and a baseline set to one's own result is the failure this package was
rebuilt to avoid (see the note on the direct-sum floor in README.md).

What is NOT obvious is whether something else overtakes that implied floor above 96, the
way Edel-Rains-Sloane overtakes Gamma_72 AT 96.  Two candidates, and both are measured here
rather than dismissed:

  the layered construction   ends at k = 24, i.e. at dimension 96, where it reaches
                             6 480 558 568 -- half the chain's total.  Beyond k = 24 there
                             is no zero-sum triple partition to run it on and no dimension
                             it could win.  cap96.py is where that is computed.

  the ERS chain              is defined for every n.  Evaluated from the same tables this
                             repository uses everywhere else, through n = 128, it exceeds
                             the implied floor by a fraction of a per cent at most.

The chain's value here is a LOWER bound on Edel-Rains-Sloane, as everywhere else in this
package -- every input is a construction.  So this script does not prove that nothing above
96 can be improved; it establishes the weaker and sufficient thing, that nothing this
repository can evaluate improves on what dimension 96 already implies, which is why 96 is
the last row rather than an arbitrary stopping point.

Exits non-zero if the excess anywhere in 97..128 reaches 1 per cent, at which point the
range would have to be reconsidered.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DELIV = os.path.join(HERE, '..', '..', '..')
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(DELIV, 'common'))

TOP = 128
THRESHOLD = 0.01          # 1 per cent: the point at which a row above 96 would be worth stating

fail = 0


def bad(msg):
    global fail
    fail += 1
    print("   FAIL: " + msg)


def floor():
    """{n: the best value implied in dimension n}, up to TOP.

    Published where a table has an entry, this repository's own claim where it does not,
    made non-decreasing, and above 96 the best direct sum of two of those.  RESULTS.md is
    the source of the claims -- imported through ers_exposure.claims(), never retyped, so
    a claim that moves moves this floor with it."""
    import ers_exposure as E
    from published import COHN
    sys.path.insert(0, DELIV)
    best = dict(COHN)                                  # 1..48 and 72, published
    for d, v in E.claims().items():                    # {dim: value}, read from RESULTS.md
        best[d] = max(best.get(d, 0), v)
    for d in range(1, TOP + 1):
        best.setdefault(d, 0)
    for d in range(2, TOP + 1):
        best[d] = max(best[d], best[d - 1])            # K is non-decreasing
    for n in range(97, TOP + 1):
        v = best[n - 1]
        for a in range(1, n // 2 + 1):
            v = max(v, best[a] + best[n - a])
        best[n] = v
    return best


def worst_excess():
    """(fraction, n): the most the chain beats the implied floor by, anywhere in 97..TOP.

    paper/factcheck.py imports this rather than parsing the table below, so the figure in
    the paper and the figure printed here have one implementation between them."""
    import ers_exposure as E
    best = floor()
    out = (0.0, None)
    for n in range(97, TOP + 1):
        exc = E.best_chain(n)[0] / float(best[n]) - 1.0
        if exc > out[0]:
            out = (exc, n)
    return out


def main():
    import ers_exposure as E

    best = floor()
    if best[96] != 12886999232:
        bad("dimension 96 is %d here, not the claimed 12886999232" % best[96])

    print("The implied floor above 96, against the Edel-Rains-Sloane chain")
    print()
    print("    n     implied floor        ERS chain      excess     the split that implies it")
    worst = (0.0, None)
    for n in range(97, TOP + 1):
        # name the split, so the floor is readable rather than a number out of a loop
        split = max(((best[a] + best[n - a], a) for a in range(1, n // 2 + 1)))
        chain = E.best_chain(n)
        if isinstance(chain, tuple):
            chain = chain[0]
        exc = chain / float(best[n]) - 1.0
        if exc > worst[0]:
            worst = (exc, n)
        if n % 4 == 1 or exc == worst[0]:
            print("  %3d %17d %16d   %+8.4f%%     %d + %d"
                  % (n, best[n], chain, 100 * exc, split[1], n - split[1]))
    print()
    if worst[1] is None:
        bad("nothing was evaluated above 96")
    elif worst[0] >= THRESHOLD:
        bad("the chain beats the implied floor by %.3f%% at n = %d, which is at or above the "
            "%.0f%% at which a row above 96 would be worth stating"
            % (100 * worst[0], worst[1], 100 * THRESHOLD))
    else:
        print("   the largest excess anywhere in 97..%d is %+.4f%%, at n = %d"
              % (TOP, 100 * worst[0], worst[1]))
        print("   so every dimension above 96 is what dimension 96 and a direct sum give,")
        print("   and the table stops at 96 because that is where it stops saying anything.")

    # ---- falsifiability: the test must be able to FAIL ---------------------------------
    # Halve the floor and the same comparison must trip.  A predicate that cannot go red is
    # a printed column, not a check.
    probe = max(E.best_chain(n)[0] / (0.5 * best[n]) - 1.0 for n in range(97, TOP + 1))
    if probe < THRESHOLD:
        bad("the comparison cannot fail: at half the floor the excess is still only %.3f%%"
            % (100 * probe))
    else:
        print("   (at half that floor the same test reads %+.1f%%, so it can fail)"
              % (100 * probe))

    print()
    print("ALL CHECKS PASSED" if not fail else "%d CHECK(S) FAILED" % fail)
    return 1 if fail else 0


if __name__ == '__main__':
    sys.exit(main())
