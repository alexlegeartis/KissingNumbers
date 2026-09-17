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
    'exhausted-18': 'the bent-coset hexagon over odd BW16 is exactly at its ceiling: tier B <= 6*M6 '
                    'and M6 = 160 (the 512-word slot graph is isomorphic to the girth-5 graph whose '
                    'independence number 160 is CP-SAT-certified), so 8358 is the scheme maximum; poles, '
                    'family size and vector shapes are all pinned (verifications/improved/dim18-bent-hexagon/README.md)',

    'exhausted-25': 'no removal is ever shared (a theorem), so the count is 196560 + H + 2 + |E| '
                    'and only the head count H moves; H = 1006 against a Delsarte ceiling of '
                    '5763, and every route beyond the template needs a code in R^24 with '
                    '~93000 points at 62 degrees (KNOWLEDGE 125, dim25-lens-heads/README)',
    'exhausted-26': 'the hexagon forces every head onto |x|^2 = 8/3 and onto the six edge midpoints, so the '
                    'layer is two head sets at cosine <= 1/4 cross-constrained at cosine 1/2; a zero-sum '
                    'triangle of norm-6 vectors gives 762 heads per side, its exact optimum, and the whole '
                    '24-cell of leans gives no more on a side in three solves, so 2 * 762 is where class '
                    'heads stop; interior heads cost >= 11 there (verify26.py; the note, section 8)',
    'exhausted-27': 'the cuboctahedron pins the twelve directions, four zero-sum triangles pairwise at '
                    '60 degrees, so the layer is four head sets at cosine <= 1/4 cross-constrained at '
                    'cosine 1/2; two coset triangles of norm-6 vectors give 2210 heads, within five '
                    'of the solver bound for that pool, and a third and fourth triangle add nothing '
                    'in an hour; the rho = 2 layer at the axis holes is empty (verify27.py; the note, '
                    'section 8)',
    'classes-28': 'the published configuration is Ma et al. 2025, arXiv:2511.13391 (PackingStar); '
                  'over it, the equator (196 560) and the direction weight (16 = 8 zero-sum '
                  'triples partitioning the 24 directions of the 24-cell, which is K(4)) are at '
                  'their ceilings, and so is the LAYER.  The weight ceiling is a THEOREM and, '
                  'here alone in 28-31, an UNCONDITIONAL one: two owner classes sharing a '
                  'direction would meet at <u,u> <= 1 everywhere and so would be one class, so '
                  'the groups are disjoint; a group is pairwise at <= -1/2, so |sum z|^2 >= 0 caps '
                  'it at three; across groups it is <= 1/2; hence the union of the groups is a '
                  '60-degree code and w = 2t + q <= floor(2 tau(k)/3) = 16, attained -- and tau(4) '
                  '= 24 is PROVED where tau(5), tau(6) and tau(7) are only best known.  Its head is a norm-8 Leech vector, not '
                  'the norm-6 of every other layer in 25-31: at height sqrt2 the cap condition '
                  'needs <v/2,u> <= sqrt6 - sqrt2, which a norm-8 head meets at <u,v> <= 2 and a '
                  'norm-6 head cannot meet at all against D4\'s covering cosine 1/sqrt2.  A head '
                  'carries every direction w with <w,w> <= 0, so the eight +-e_k of R^4; two heads '
                  'share a direction iff <v,v> <= 0, so a frame and its negatives share all eight; '
                  'and 48 pairwise non-positive vectors in R^24 ARE a frame, so H <= 24 head lines '
                  'and 16H = 384 points is the maximum.  The axis falls from 24 to the 16 '
                  'half-vectors of the dual 24-cell, so the layer nets 16H - 8 = +376.  That axis '
                  'is MAXIMAL, exhaustively: a 17th point is a UNIT vector of the polytope cut out '
                  'by <a,a_i> <= 1/2, <a,z> <= sqrt3/2 and |<a,w>| <= 1/sqrt2, every right-hand '
                  'side is positive so the origin is interior, |a| is convex so its maximum is at '
                  'a VERTEX, and all 48 vertices have |a| <= 0.765 -- no unit vector exists at all '
                  '(ceilings.py; the control is that dropping any one of the 16 recovers it, and '
                  'nothing else).  What is left is only the class: +32 per extra LINE (weight 16, two points per weight '
                  'unit per line) against 248 realised and 425 by the Delsarte LP, and here the '
                  'class must additionally consist of type-B lines of shape (2^8, 0^16) in the '
                  'frame coordinates, which is what makes every owner meet every head at |<u,v>| '
                  '<= 2.  Such classes exist: the six tetrads of a sextet carry 24 mutually '
                  'orthogonal sign-lines forming one Leech frame, and in it the record 248-line '
                  'class is all type-B (verifications/improved/dim28-norm8-frame-layer/README.md)',
    'classes-29': 'the published configuration is Ma et al. 2025, arXiv:2511.13391 (PackingStar).  '
                  'Over it the equator (196 560), the direction weight (26 = 12 zero-sum triangles '
                  'plus 2 antipodal pairs) and the LAYER are at their ceilings.  The weight ceiling '
                  'is a THEOREM -- disjoint groups, at most three directions each, union a '
                  '60-degree code, so w = 2t + q <= floor(2 tau(k)/3) = 26 -- but a CONDITIONAL '
                  'one, since tau(5) = 40 is the best known and not a proved value: against the '
                  'proved tau(5) <= 44 the unconditional ceiling is 29, worth 1 488 points, and '
                  'reaching it needs a better 5-dimensional kissing number.  The layer is 480 = 48 frame '
                  'vectors x 10 cross-polytope directions, which is every head line of a Leech frame '
                  'on every direction a head can carry; nothing at height sqrt2 can be larger.  The '
                  'packing is exact: 14 pairwise disjoint classes of the full 248 lines.  Two things '
                  'move it.  (a) The class: +52 per extra LINE (weight 26, two points per weight unit '
                  'per line), against 248 realised and 425.45 by the Delsarte LP on the Leech line '
                  'scheme -- a Leech question, shared with dimensions 28, 30, 31, 49-61 and 73-95.  '
                  '(b) tau(5) itself, above.  The AXIS is not a lever: 32 of the published 40 '
                  'survive |<a,w>| <= 1/sqrt2 -- the 8 that die are at |<a,w>| = 1 EXACTLY, they lie '
                  'along a frame direction -- and those 32 are MAXIMAL exhaustively, not by '
                  'sampling: the admissible set is a polytope with the origin interior, |a| is '
                  'convex so its maximum is at a VERTEX, and all 114 vertices have |a| <= 0.7929, so '
                  'no unit vector can be added at all and no one-for-one swap opens it either '
                  '(ceilings.py 29).  Only a different 60-degree code of R^5 would change that '
                  '(research/collab2531/FINDINGS.md section 59; dim29-30-frame-layer/README.md)',
    'classes-30': 'the equator (196 560), the direction weight (48 = floor(2 x 72 / 3)) and the '
                  'LAYER are at their ceilings.  The weight ceiling is a THEOREM -- disjoint '
                  'groups, at most three directions each, union a 60-degree code, so '
                  'w = 2t + q <= floor(2 tau(k)/3) -- but a CONDITIONAL one, since tau(6) = 72 is '
                  'the best known and not a proved value: against the proved tau(6) <= 77 the '
                  'unconditional ceiling is 51, worth 1 488 points.  The layer is '
                  '576 = 48 frame vectors x 12 cross-polytope directions, every head line of a Leech '
                  'frame on every direction.  Two things move it.  (a) The class: +96 per extra LINE, '
                  'against 248 realised and 425.45 by the Delsarte LP -- the same Leech question as in '
                  'dimensions 28, 29, 31, 49-61 and 73-95.  The PACKING is CLOSED: the 24 classes are '
                  'pairwise disjoint and full, 24 x 248 = 5 952 owner lines, found by freezing a core of 12 '
                  'disjoint classes, building a pool of monomial images disjoint from the whole core, and '
                  'solving max-clique on that pool -- the coordinate descent it replaces was stuck at 5 941 '
                  '(research/collab2531/BRAINSTORM_2831.md).  (b) The axis is now 68 of 72, not the 48 the published axis kept.  A maximum '
                  '60-degree code in R^6 is a ROTATED E6 and the rotation is free: 45, 45 and 90 '
                  'degrees in the coordinate planes (0,1), (2,3), (3,5) keeps every coordinate in the '
                  '(a + b r2 + c r3 + d r6)/24 field the package already stores, and 68 of the 72 are '
                  'flat enough for the layer.  That 68 is MAXIMAL, and exhaustively so rather than '
                  'by the 6e6 samples first used: the admissible set is a polytope with the origin '
                  'interior, |a| is convex so its maximum is at a VERTEX, and all 202 vertices have '
                  '|a| <= 0.8165, so no unit vector can be added at all and no one-for-one swap '
                  'opens it either (ceilings.py 30).  The last 4 would need a different axis code '
                  'altogether.  So the levers are the class and tau(6) '
                  '(research/collab2531/BRAINSTORM_2831.md; dim29-30-frame-layer/README.md)',
    'classes-31': 'the equator (196 560), the direction weight (84 = floor(2 x 126 / 3)), the AXIS '
                  'and the LAYER are at their ceilings.  The weight ceiling is a THEOREM -- '
                  'disjoint groups, at most three directions each, union a 60-degree code, so '
                  'w = 2t + q <= floor(2 tau(k)/3) -- but a CONDITIONAL one, since tau(7) = 126 is '
                  'the best known and not a proved value: against the proved tau(7) <= 134 the '
                  'unconditional ceiling is 89, worth 2 480 points.  The axis of 118 is MAXIMAL '
                  'exhaustively: the admissible set is a polytope with the origin interior, |a| is '
                  'convex so its maximum is at a VERTEX, and all 1120 vertices have |a| <= 0.8660, '
                  'so no unit vector can be added at all (ceilings.py 31).  The layer is 672 = 48 frame vectors x '
                  '14 cross-polytope directions, every head line of a Leech frame on every direction, '
                  'and it deletes no equator point.  Two things move it.  (a) The class: +168 per extra '
                  'LINE, the largest per-line value of any dimension, against 248 realised and 425.45 '
                  'by the Delsarte LP -- the same Leech question as in dimensions 28, 29, 30, 49-61 and '
                  '73-95.  (b) The PACKING: the 42 classes carry 10 328 of the 10 416 lines that 42 full '
                  'classes would give, so closing the remaining 88 repeats is worth up to +352.  That '
                  'is a search, not an obstruction: two random monomial images are disjoint with '
                  'probability 0.4319 +- 0.0030 (0.4155 was measured with a biased sampler; see '
                  '145; re-measured independently over 40 000 uniform images at 0.4305 +- 0.0025), '
                  'so the random-graph clique number on the ~1e12 images is '
                  '2 ln(1e12)/ln(1/0.4319) = 66 and 42 is far below it.  What is closed is the shortcut -- no monomial element '
                  'acts on the type-B lines with order >= 42, and C_M24(7A) = C_7 x S_3, the order-42 '
                  'group that does exist, leaves only 92.7 per cent of the lines in a free orbit, so a '
                  '248-line transversal has probability 7e-9 x 7.6e-13.  Two NEIGHBOURHOODS are '
                  'now closed exhaustively rather than by a search stopping: the sign words, '
                  'where every one of the 861 exact pair tables reaches 0 somewhere but annealing '
                  'cannot leave 88; and the RE-SLOTTING, where one trio admits exactly 2016 '
                  'supports and 64512 x 8 / 2016 = 256 images realise each, so replacing every '
                  'class by each of its 256 same-support images over all 4096 sign words is 44 '
                  'million placements and zero improvements.  The collisions are also completely '
                  'named: every slot lies in a coset of one of the 30 extended Hamming codes of '
                  'its octad, a sign word moves the COSET and never the CODE, and the code is a '
                  'function of the class\'s TRIO -- it IS one of the 15 trios through that octad '
                  '-- so two classes there collide with probability 1/8 when the trios agree and '
                  'exactly 1/2 when they do not, and the 120 code pairs meeting in dimension 0, '
                  'where no sign word could separate them, never occur.  None of that closes the '
                  '88: on octad coverage, bit budget, pair probability and code agreement alike, '
                  'five fresh builds match the shipped family and lose twice as many lines, so '
                  'what is left is search depth '
                  '(research/collab2531/BRAINSTORM_2831.md section 9; '
                  'research/collab2531/dim31/{kappa31,reslot31}.py; '
                  'dim31-frame-layer/README.md)',
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
for _d, _s, _r in ((18, 'exhausted', 'exhausted-18'),
                   (25, 'exhausted', 'exhausted-25'),
                   (26, 'exhausted', 'exhausted-26'),
                   (27, 'exhausted', 'exhausted-27'),
                   (28, 'classes', 'classes-28'),
                   (29, 'classes', 'classes-29'),
                   (30, 'classes', 'classes-30'),
                   (31, 'classes', 'classes-31'),
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
