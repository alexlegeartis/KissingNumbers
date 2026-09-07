#!/usr/bin/env python3
"""What would have to change for each claim to improve.  The status column of RESULTS.md.

    from status import STATUS, LEGEND, status_for

READ THIS FIRST.  A status describes THIS PROJECT'S MECHANISM in that dimension.  None of them
says anything about tau(n) itself.  No matching upper bound is known in any dimension this
repository touches, so no claim here is known to be optimal -- "exhausted" means the idea is
finished, not the dimension.

The point of the column is that 47 dimensions collapse into four actions:

  exhausted   The mechanism provably cannot give more here.  Do not retry the idea; a
              different construction is needed.
  external    The value is a function of a PUBLISHED TABLE.  It improves the moment that table
              does, with no new mathematics -- so recheck it when the table moves.  Dimensions
              62 and 63 are here because the Edel-Rains-Sloane chain beats the cap
              construction's CEILING there, not merely the family it happens to have.
  classes     The direction side and the axis side are closed and only the class family is
              short of its ceiling.  Every dimension with this status improves TOGETHER when a
              better family of pairwise-disjoint classes is found, which is where the compute
              goes.
  lp-exact    The value is the EXACT OPTIMUM of the k-point programme for a k-tuple Gram
              whose determinant is at the Hermite floor: a rank-k sublattice of Gamma_72 has
              det >= (8/gamma_k)^k, with equality only for the densest rank-k lattice, and
              for k = 1..4 the attaining Gram (A_2, A_3, D_4) is exhibited.  The floor is a
              theorem.  Optimum, not merely a certificate: common/exact_vertex.py solves the
              dual vertex in rational arithmetic and verifies it over every column, and the
              numeric primal attains the same value, so nothing about this programme improves
              it.  That the lowest determinant gives the best bound is a MEASUREMENT across
              the Grams tried; the basis is not a variable at all, all sixteen k = 3
              presentations giving one answer.

The reasons are one line each and cite the script or section that establishes them.
"""
import os

from published import TAU

LEGEND = {
    'exhausted': 'the mechanism cannot give more; a different construction is needed',
    'external': 'a function of a published table -- recheck when the table moves',
    'classes': 'only a better class family helps; these all improve together',
    'lp-exact': 'the programme solved to its EXACT optimum, for a Gram at the Hermite determinant floor; only a larger k helps',
}

# The largest class actually shipped, and what it is worth against Caro-Wei.  This used to
# be the literal 2111; the greedy family was rebuilt on 2026-08-26 and the largest class
# became 2118, so the sentence described a class the repository no longer held.  It is read
# from the shipped size tables instead, and falls back to the values as of that date when the
# data directory is not present (a clean checkout of the paper alone).
_CW_G72, _CW_LEECH = 221, 11.3


def _g72_best():
    import array
    best = 0
    d = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'verifications',
                     'improved', 'dim73-95-gamma72-caps', 'data')
    for fn in ('class_sizes_32000.npy', 'class_sizes_gpu.npy'):
        f = os.path.join(d, fn)
        if not os.path.exists(f):
            continue
        import ast
        with open(f, 'rb') as fh:
            assert fh.read(6)[1:] == b'NUMPY'   # the magic, without the 0x93 byte
            major = ord(fh.read(1))
            fh.read(1)
            n = int.from_bytes(fh.read(2 if major == 1 else 4), 'little')
            h = ast.literal_eval(fh.read(n).decode('latin1').strip())
            a = array.array({'|i1': 'b', '<i2': 'h', '<i4': 'i', '<i8': 'q',
                             '<u2': 'H', '<u4': 'I', '<u8': 'Q'}[h['descr']])
            a.frombytes(fh.read())
        best = max(best, max(a))
    return best or 2118


def _axis_layers():
    """{k: poles} for the shipped axis layers, from the .npz headers -- stdlib only.

    Quoted here rather than computed, these numbers went stale twice in one day."""
    import zipfile
    out = {}
    d = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                     'verifications', 'improved', 'dim73-95-gamma72-caps', 'data')
    for k in range(1, 24):
        # FIVE files can hold a layer at one k: from the record configuration of R^k,
        # from the same idea over Z[sqrt2], from a rotated copy of a Kappa section, and from
        # a rotated copy of the cap directions under either naming.  The largest is the one
        # the dimension ships, which is how final81-95.py chooses; reading only the rotated
        # file under-reports three of them, and the Kappa pattern was missing until
        # 2026-08-30, which hid k = 12 and 13 entirely.
        for nm in ('poles_rec_%d.npz' % k, 'poles_sqrt2_%d.npz' % k,
                   'poles_kap_%d.npz' % k,
                   'poles_rot_%d.npz' % k, 'poles_rot_k%d.npz' % k):
            f = os.path.join(d, nm)
            if not os.path.exists(f):
                continue
            try:
                with zipfile.ZipFile(f) as z:
                    with z.open('keep.npy') as h:
                        head = h.read(128)
                i = head.index(b'(') + 1
                out[k] = max(out.get(k, 0), int(head[i:head.index(b',', i)]))
            except Exception:
                pass
    return out


def _axis_facts():
    """(best layer, its k, the shell layer there, the k at which the layer reaches tau,
    the total shortfall below tau over every k that has a layer)

    The shortfall is what is left in the axis layers, and it is the figure the sentence
    below used to give as "a few tens of points".  It is 602, which is a few hundred -- so
    it is summed here rather than described."""
    lay = _axis_layers()
    if not lay:
        return 0, 0, 0, '', 0
    kbest = max(lay, key=lambda k: lay[k])
    shell = 0
    d = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                     'verifications', 'improved', 'dim73-95-gamma72-caps', 'data')
    for nm in ('poles_%d.npy' % kbest, 'poles_k%d.npy' % kbest):
        f = os.path.join(d, nm)
        if os.path.exists(f):
            with open(f, 'rb') as fh:
                import ast
                assert fh.read(6)[1:] == b'NUMPY'
                major = ord(fh.read(1))
                fh.read(1)
                n = int.from_bytes(fh.read(2 if major == 1 else 4), 'little')
                shell = 2 * ast.literal_eval(fh.read(n).decode('latin1').strip())['shape'][0]
    at = sorted(k for k in lay if k in TAU and lay[k] == TAU[k])
    short = sum(TAU[k] - lay[k] for k in lay if k in TAU and lay[k] < TAU[k])
    return (lay[kbest], kbest, shell, ', '.join('k = %d' % k for k in at) or 'no k', short)


_B = _g72_best()
_AX = _axis_facts()
_G72_BEST = (_B, _B / float(_CW_G72), 100.0 * (_CW_LEECH * _CW_G72 / float(_B) - 1.0))

REASON = {
    'exhausted-25': 'no removal is ever shared (a theorem), so the count is 196560 + H + 2 + |E| '
                    'and only the head count H moves; H = 1006 against a Delsarte ceiling of '
                    '5763, and every route beyond the template needs a code in R^24 with '
                    '~93000 points at 62 degrees (KNOWLEDGE 125, dim25-lens-heads/README)',
    'exhausted-27': 'four classes suffice where the published configuration used five, and the '
                    'configuration is verified over all 20 108 045 530 pairs '
                    '(scripts/verify_exhaustive.py)',
    'exhausted-38': 'BOTH layers are at their ceilings: the cap layer at 3 x 196 560 = '
                    '589 680 with the equator empty, and the axis layer at a full tau(14) = '
                    '1932, so the total 591 612 IS 3 tau(24) + tau(14) and this mechanism has '
                    'no slack left at all.  It read 1908 until 2026-08-27, when the axis '
                    'rotation was found by DESCENT on SO(14) rather than by sampling -- '
                    'sampling loses 56 on average and cannot get there, which is a '
                    'measurement, but "so it is not there" was not.  The one caveat is that '
                    'tau(14) = 1932 is the best KNOWN dimension-14 kissing number, not a '
                    'proven one, so improving dimension 14 improves this by the same amount '
                    'and nothing else does.  The construction was run in every dimension '
                    '32-40 and wins only here (dim38/README, scripts/axis_rotation.py)',
    'external-6263': 'the Edel-Rains-Sloane chain (n,15,2): level 0 is the sign code [62,26,16] '
                     'or [63,27,16] from codetables.de, levels 1 and 2 are built here and cite '
                     'nothing; the P_48 cap construction cannot reach even level 0 with ANY '
                     'class family -- its ceiling is 66075240 and 70543480 '
                     '(dim49-63-p48-caps/scripts/ceiling.py)',
    'external-96': 'the Edel-Rains-Sloane chain (96,24,6,1).  Level 0 is the sign code '
                   '[96,33,24] from codetables.de, and its 2^33 is two thirds of the total; '
                   'the support codes are the grid maps [24,9,12]_4 and RS[6,4,3]_16, the '
                   'latter built here.  It improves whenever any of those tables does.  The '
                   'Gamma_72 cap construction is not the rival: at k = 24 it reaches only '
                   '6 480 558 568, so 96 is out of its reach entirely, which is why this '
                   'package sits in closed/ (dim96-ers-takeover/)',
    'external-39': 'the Edel-Rains-Sloane chain evaluated with current constant-weight tables; '
                   'the whole +128 came from A(39,8,8) moving 3323 -> 3324, so it improves '
                   'whenever Brouwer\'s table does',
    'classes-p48': 'directions and axis layer are at their ceilings, and the class family is '
                   'close to its own: scripts/ceiling.py shows that 94-105% of the distance '
                   'from the family to the ceiling is the overlap two random automorphism '
                   'images have BY CHANCE, so best-of-N selection has almost nothing left to '
                   'win -- a real gain needs a family disjoint BY CONSTRUCTION (KNOWLEDGE 89, '
                   '91).  The same script shows the ceiling is BELOW the level-0 sign code in '
                   'dimensions 62 and 63, which is why they are not in this group',
    'classes-g72': ('the cap DIRECTIONS are at their ceiling -- the partition into zero-sum '
                   'triples is perfect where one exists, and three is the most a class line '
                   'can ever serve.  The class FAMILY is what moves, and it moved twice in '
                   'one session (KNOWLEDGE 90, 91): the class size %d is a factor %.1f over '
                   "Caro-Wei against the Leech lattice's 11.3, so matching that would be a "
                   'further %.0f%% on every one of these 23 dimensions.  The AXIS layer was '
                   'rebuilt on 2026-08-27 as a ROTATED COPY of the cap directions, which is a '
                   '60-degree code for free and only has to clear the 30-degree caps: at '
                   'k = %d it gives %d poles against the shell rule\'s %d, and it reaches '
                   'tau(k) exactly at %s.  What is left '
                   'there is the handful a rotation loses to the caps, and, where the cap set '
                   'has fewer than tau(k) vectors, the gap between the two -- %d points in '
                   'total, against a claim of 6.4e9, and 321 of them are the two the record '
                   'configuration of R^k cannot be carried into (k = 12 is given in floating '
                   'point by the data set, k = 13 differs from K_13 in Hasse invariant).  The '
                   'class family is the only thing that moves this group '
                   '(scripts/axis_rotate.py, scripts/axis_block.py, '
                   'scripts/verify_poles.py)') % (_G72_BEST + (_AX[1], _AX[0], _AX[2],
                                                               _AX[3], _AX[4])),
    'lp-exact-g72': 'the k-point programme solved to its EXACT optimum by '
                    'common/exact_vertex.py, which substitutes out the identifications and '
                    'solves the dual vertex in rational arithmetic rather than rounding a '
                    'float, so the value is the most this programme can certify: no basis, '
                    'dual or solver setting moves it, and all sixteen determinant-256 '
                    'presentations agree at k = 3.  The Gram determinant is at the Hermite '
                    'floor, which is a theorem, and the tuple is exhibited in Gamma_72.  That '
                    'the lowest determinant gives the best bound is measured across the k = 3 '
                    'table, not proved.  So only a larger k could help, and k = 5 extrapolates '
                    'below the floor (scripts/verify68.py, verify69.py, k4_gram_min_det.py, '
                    'presentation_invariance.py)',
}

STATUS = {}
for _d, _s, _r in ((25, 'exhausted', 'exhausted-25'),
                   (27, 'exhausted', 'exhausted-27'),
                   (38, 'exhausted', 'exhausted-38'),
                   (39, 'external', 'external-39')):
    STATUS[_d] = (_s, REASON[_r])
for _d in range(49, 62):
    STATUS[_d] = ('classes', REASON['classes-p48'])
for _d in (62, 63):
    STATUS[_d] = ('external', REASON['external-6263'])
for _d in (68, 69, 70, 71):
    STATUS[_d] = ('lp-exact', REASON['lp-exact-g72'])
for _d in range(73, 96):
    STATUS[_d] = ('classes', REASON['classes-g72'])
STATUS[96] = ('external', REASON['external-96'])


def runs(ds):
    """compress a dimension list into contiguous runs, e.g. [49..61, 73..95] -> "49-61, 73-95"

    'min-max' is wrong as soon as a status is not an interval, and the class statuses are not:
    they are 49-61 and 73-95, and printing "49-95" silently claims 62-72 as well."""
    out = []
    for d in sorted(ds):
        if out and d == out[-1][1] + 1:
            out[-1][1] = d
        else:
            out.append([d, d])
    return ', '.join('%d-%d' % (a, b) if b > a else str(a) for a, b in out)


def status_for(d):
    """(status, reason) for a dimension, or (None, None) if it has not been classified"""
    return STATUS.get(d, (None, None))


if __name__ == '__main__':
    print(__doc__)
    by = {}
    for d, (s, _r) in sorted(STATUS.items()):
        by.setdefault(s, []).append(d)
    for s in sorted(by):
        ds = by[s]
        print('  %-10s %2d dimensions  %-16s %s' % (s, len(ds), runs(ds), LEGEND[s]))
