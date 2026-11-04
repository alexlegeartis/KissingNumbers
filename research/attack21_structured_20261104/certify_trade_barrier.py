"""Solver-independent certificate: no gain removing at most twelve old holes.

This applies to ANY new dominating set, including ones with adjacent holes.
Thus the local exclusion covers the full line-label model, not just the
independent-domination sufficient construction.
"""
import json
import numpy as np
from model import C, N, SIDE, DATA


def parameters(H):
    count = H[N].sum(1)
    assert H.sum() == 64 and not count[H].any()
    private = (~H) & (count == 1)
    assert np.all(private[N[H]].sum(1) == 12)
    assert count[~H].max() == 3
    # Every new hole covers at most seventeen private vertices by adjacency.
    assert private[N[~H]].sum(1).max() == 17
    # Its closed neighborhood can cover at most two private neighbors of any
    # single removed old hole. For same-side vertices closure adds nothing.
    owners = N[private][H[N[private]]]
    old = np.flatnonzero(H)
    owner_index = {int(u):j for j,u in enumerate(old)}
    P = np.zeros((1024,64),dtype=np.int16)
    for u,v in zip(np.flatnonzero(private),owners): P[u,owner_index[int(v)]]=1
    common = P[N].sum(1)
    assert common[~H].max() == 2
    return dict(holes=64,private_per_hole=12,max_old_neighbors=3,
                max_private_neighbors=17,max_private_from_one_hole=2)


def admissible(a,b,c,d):
    # a,b removed old holes; c,d added new holes, by bipartition.
    # Removed A holes have 12a private vertices. New A holes cover at most
    # min(2ac,17c) of these through edges; new B holes cover at most d via self.
    return (12*a <= min(2*a*c,17*c)+d and
            12*b <= min(2*b*d,17*d)+c and a <= 3*d and b <= 3*c)


def tuples(radius):
    return [(a,r-a,c,d) for r in range(1,radius+1) for a in range(r+1)
            for c in range(r) for d in range(r-c) if admissible(a,r-a,c,d)]


if __name__ == '__main__':
    seeds = np.concatenate([np.load(DATA/'recovery.npz')['holes'],
                            np.load(DATA/'nonlinear_quotient_holes.npy')])
    checks = [parameters(H) for H in seeds]
    possible = tuples(13)
    assert not tuples(12)
    assert possible and all(a+b==13 for a,b,c,d in possible)
    result = dict(status='PASS', kernels_checked=len(seeds), removed_holes_excluded=12,
                  added_holes_maximum='one fewer than removed',
                  scope='Any smaller dominating hole set; hence all line labelings near these seeds.',
                  arithmetic='exhaustive small integer inequalities; no numerical or SAT tolerance',
                  parameters=checks[0], first_unexcluded_tuples=possible)
    (DATA/'trade_barrier_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
