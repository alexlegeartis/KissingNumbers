#!/usr/bin/env python3
"""The published state of the art, with provenance.  Single source of truth for the
"previously known" column of every claim in this repository.

    from published import COHN, FLOOR, floor_for, REFERENCES

COHN
    The lower-bound column of H. Cohn, *Table of kissing number bounds*,
    <https://cohn.mit.edu/kissing-numbers>, permanent copy
    <https://hdl.handle.net/1721.1/153312>.  Re-fetched 2026-08-20, 2026-08-25, again
    2026-08-26 when dimension 96 was claimed, and again 2026-08-30 on the day of release,
    when dimensions 25, 27, 38, 39, 44, 45, 46, 47, 48 and 72 -- every row any claim or
    withdrawal here depends on -- were read off the live page and all ten agreed.
    On 2026-08-26 ALL 49 entries were compared,
    not a sample, and every one agreed; the upper-bound column was compared too and its 43
    genuine entries agreed (the six rows that appear to disagree are exactly the solved
    dimensions 1, 2, 3, 4, 8 and 24, where the ratio column reads 1).  So the four claims
    that fall inside the table (25, 27, 38, 39) are measured against current values.
    THE TABLE MOVES: two entries (dimensions 12 and 19) changed while this project was in
    progress.  Re-check before claiming anything.

    The table covers dimensions 1-48 and 72, AND NOTHING ELSE.  This matters enormously
    above dimension 48 and is the single fact most often got wrong.

FLOOR
    Above dimension 48 there is no table, so "the value to beat" is what monotonicity of
    tau inherits from the nearest tabulated dimension below, or a direct sum, whichever is
    larger.  floor_for(d) returns (value, explanation).  The three regimes:

      49-63   52416000 + tau(d-48)   direct sum P_48 (+) (best (d-48)-dimensional lattice
                                     rescaled to minimum 6).  Cohn's table has no entry in
                                     this range and neither did Nebe-Sloane, so this is the
                                     floor -- but it is NOT the whole story, and the way it
                                     was got wrong twice is worth keeping.

                                     The argument used to be "the value 2^28 = A(64,16) that
                                     carries dimension 64 is not available below n_0 = 64".
                                     The premise is true and the conclusion does not follow:
                                     2^27 is available at n_0 = 63 and 2^26 at n_0 = 62, from
                                     [63,27,16] and [62,26,16] with ceil(63/4) = ceil(62/4)
                                     = 16, so a SINGLE SIGN CODE beats the cap construction
                                     there.  That was found on 2026-08-23 and the two claims
                                     moved to verifications/improved/dim62-63-ers-chain/.
                                     They could not be rescued: the cap construction's
                                     absolute ceiling, every one of its lambda(k) classes at
                                     the maximum size 7069, is 66075240 and 70543480, both
                                     short (dim49-63-p48-caps/scripts/ceiling.py).

                                     Then the SAME mistake one level up.  The sign code is
                                     only LEVEL 0 of Edel-Rains-Sloane.  The full chain adds
                                     A(n, n_1, n_1) * A(n_1, ceil(n_1/4)) and more, which is
                                     what carries their own dimension-64 total from 2^28 =
                                     268435456 to 331737984.  Evaluated at 62 and 63 with the
                                     chain (n, 15, 2) it gives 71310732 and 138419844, and
                                     those are what is claimed now.  Levels 1 and 2 are BUILT
                                     in that package, not cited.  Found 2026-08-25 by
                                     verifications/closed/dim96-ers-takeover/ers_exposure.py,
                                     which evaluates the full chain in EVERY claimed dimension
                                     rather than only level 0 as ers_sweep.py does.

                                     Dimensions 49-61 are clear: the full chain reaches
                                     18400096 at n = 49-59 (chain (48,8,2), the best sign code
                                     being A(48,12) = 2^24) and 37755816 at n = 60-61 (chain
                                     (60,15,2), A(60,15) = 2^25), against claims of 52.4M to
                                     59.9M.  The authors themselves say that at dimension 48
                                     their construction "gives less than half" of 52416000.
      64-71   331737984              dimension 64, Edel-Rains-Sloane 1998, inherited by
                                     monotonicity.  Their construction is generic in n but was
                                     evaluated by its authors only at 32, 36, 40, 44, 64, 80,
                                     128.  Evaluated at 65-71 with the chain (64,16,4,1) and
                                     the only citable A(n,16,16) -- their own 30828, monotone
                                     in n -- it gives slightly MORE than 331737984, because
                                     A(n,4,4) and A(n,1,1) grow: 331743186 at n = 65 up to
                                     331796542 at n = 71 (A(68,4,4) = C(68,3)/4 exactly, since
                                     68 = 2 mod 6 admits an S(3,4,68)).  Every claim here
                                     clears that, dimension 68 by 8.9% and the rest by 1.9x or
                                     more.  Generously handing them the LENGTH-80 constant-
                                     weight code at length 70 or 71 gives at most 563125646,
                                     which dimensions 70 and 71 also clear and dimension 68
                                     does NOT -- so dimension 68's margin rests on nobody
                                     having published an A(68,16,16), and it would fall to one
                                     of size 45235.  That is the most exposed claim in the
                                     repository and ers_exposure.py ranks it as such.
      73-95   6218175600 + tau(d-72) direct sum Gamma_72 (+) (best (d-72)-dimensional
                                     lattice rescaled to minimum 8).  This is exactly the
                                     kind of entry the old Nebe-Sloane table recorded at
                                     dimension 80 as "Gamma_72 perp E8".  Edel-Rains-Sloane
                                     stays far below: its dimension-80 total is 1368532064,
                                     and codetables.de gives at most A(n_0, n_0/4) = 2^30
                                     all the way to n_0 = 92, and 2^31, 2^31, 2^32 at
                                     n_0 = 93, 94, 95.  A(96,24) >= 2^33 ([96,33,24]) is
                                     where that stops: ERS at n = 96 totals 12886999232,
                                     a factor 2.07 on Gamma_72, which is why this
                                     repository ends at 95.  Worked out in
                                     verifications/closed/dim96-ers-takeover/.

    The Nebe-Sloane "Table of the highest kissing numbers presently known" that Cohn's
    table replaced (last modified Feb 2011, retrievable from the Wayback Machine) listed
    dimensions 1-40, 42, 44, 48, 64, 72, 80, 128 and nothing else, so it does not fill the
    gaps either.  OEIS is no help: A001116 stops at n = 9, A002336 at 25, A028923 at 17.

    python published.py    prints the table with its provenance.
"""

# ---------------------------------------------------------------------------------------
# H. Cohn, Table of kissing number bounds, lower-bound column.  Verified 2026-08-20,
# re-fetched 2026-08-25, 2026-08-26 and 2026-08-30, unchanged at every row this repository
# reads.
# The bracketed keys are the reference labels of REFERENCES below.
COHN = {
    1: 2, 2: 6, 3: 12, 4: 24, 5: 40, 6: 72, 7: 126, 8: 240, 9: 306, 10: 510, 11: 604,
    12: 841, 13: 1154, 14: 1932, 15: 2564, 16: 4320, 17: 5730, 18: 7654, 19: 11948,
    20: 19448, 21: 29768, 22: 49896, 23: 93150, 24: 196560, 25: 197056, 26: 198550,
    27: 200044, 28: 204520, 29: 209496, 30: 220440, 31: 238350, 32: 345408, 33: 360640,
    34: 380868, 35: 409548, 36: 484568, 37: 494312, 38: 566652, 39: 755988, 40: 1064368,
    41: 1170384, 42: 1250676, 43: 2060399, 44: 2948552, 45: 3047160, 46: 5318060,
    47: 9741412, 48: 52416000, 72: 6218175600,
}

# The Nebe-Sloane "Table of the highest kissing numbers presently known" that Cohn's
# table replaced.  Kept as data and not only as prose because ONE of these dimensions,
# 80, lies inside a range this repository claims in, and both the paper and factcheck
# had asserted that no table has an entry anywhere above 48.  The entry there is
# Gamma_72 perp E_8 and it agrees with floor_for's direct sum to the digit -- asserted
# below, once floor_for exists.  Anything that wants to say "untabulated" must test
# against BOTH tables, not against COHN alone.
NEBE_SLOANE = tuple(list(range(1, 41)) + [42, 44, 48, 64, 72, 80, 128])

# who holds each of the entries this repository touches or cites
COHN_SOURCE = {
    11: "Bianchi-Kwon-Pappu-Zou 2026, arXiv:2606.10402",
    12: "Takhanov-Assylbekov-Yun 2026, arXiv:2606.18984",
    13: "Zinoviev-Ericson 1999",
    17: "Cohn-Li 2024, arXiv:2411.04916",
    18: "Cohn-Li 2024, arXiv:2411.04916",
    19: "Ho 2026, arXiv:2603.10425",
    20: "Cohn-Li 2024, arXiv:2411.04916",
    21: "Cohn-Li 2024, arXiv:2411.04916",
    22: "Leech-Sloane 1971 (Lambda_22)",
    23: "Leech-Sloane 1971 (Lambda_23)",
    24: "Leech 1967 (the Leech lattice); optimality Odlyzko-Sloane 1979, Levenshtein 1979",
    25: "Ma et al. 2025, arXiv:2511.13391 (PackingStar)",
    26: "Ma et al. 2025, arXiv:2511.13391",
    27: "Ma et al. 2025, arXiv:2511.13391",
    28: "Ma et al. 2025, arXiv:2511.13391",
    29: "Ma et al. 2025, arXiv:2511.13391",
    30: "Ma et al. 2025, arXiv:2511.13391",
    31: "Ma et al. 2025, arXiv:2511.13391",
    38: "Edel-Rains-Sloane 1998 with Brouwer's constant-weight code table",
    39: "Edel-Rains-Sloane 1998 with Brouwer's constant-weight code table",
    43: "Sun-Wang 2026, arXiv:2607.20359 (antipode construction)",
    44: "Edel-Rains-Sloane 1998",
    45: "Edel-Rains-Sloane 1998 with Brouwer's constant-weight code table",
    46: "Edel-Rains-Sloane 1998 with Brouwer's constant-weight code table",
    47: "Edel-Rains-Sloane 1998 with Brouwer's constant-weight code table",
    48: "Leech-Sloane 1971 (P_48)",
    72: "Nebe 2012, J. Reine Angew. Math. 673, 237-247 (Gamma_72)",
}

# The best known kissing numbers used as POLE configurations by the cap construction, i.e.
# tau(k) for small k.  Same table as COHN, repeated here so the cap packages need only one
# import.  Entries 17-21 are Cohn-Li, 19 is Ho, 12 is Takhanov et al., 11 is Bianchi et al.
TAU = {k: COHN[k] for k in range(1, 24)}

# Values published outside Cohn's table that act as ceilings or floors above dimension 48.
# Edel-Rains-Sloane 1998, section 3: their seven worked examples, and there are no others.
# The paper is titled "in dimensions 32 to 128", but that is the SPAN of these seven, not a
# table by dimension -- it tabulates no values and does not mention 96.  All seven reproduce
# from their own stated inputs; closed/dim32-44-ers-audit/ checks 64 and 80 against the
# formula.  Checked against arXiv:math/0207291 on 2026-08-26.
#
# Kept here so that no floor can fall below a published value.  DIMENSION 80 IS THE TRAP: it
# lies inside 73-95, and three files in this repository asserted that nobody had published ERS
# anywhere in that range.  The conclusion survived -- the direct sum 6218175840 beats their
# 1368532064 by 4.5x -- but the reason given for it was false.
ERS_PUBLISHED = {32: 276032, 36: 438872, 40: 991792, 44: 2948552,
                 64: 331737984, 80: 1368532064, 128: 8863556495104}

# Values above dimension 48 that are not in any table.  READ THE COMMENT ON EACH: these four
# entries are four DIFFERENT KINDS of number, and two of them are traps.  Nothing may treat
# this dict as "the best known value in dimension d" -- only floor_for() decides that.
OTHER = {
    64: 331737984,          # BEST KNOWN.  Edel-Rains-Sloane 1998, non-lattice; one of their
                            # seven worked examples, and floor_for uses it for 64-71.
    80: 6218175600 + 240,   # BEST KNOWN.  Nebe-Sloane table, "Gamma_72 perp E8"; the same
                            # value floor_for's 73-95 branch computes as a direct sum.  ERS
                            # also published dimension 80, at 1368532064, and this beats it
                            # by 4.5x -- see ERS_PUBLISHED above.
    96: 12886999232,        # NOT PUBLISHED, AND NOT A FLOOR.  This is THIS REPOSITORY'S OWN
                            # CLAIM in dimension 96, the Edel-Rains-Sloane chain (96,24,6,1)
                            # evaluated with current code tables.  It was floor_for(96) until
                            # 2026-08-26, which hid the claim by measuring it against itself.
                            # Kept only for the monotonicity message below.  The floor at 96
                            # is the direct sum; see floor_for's docstring.
    128: 218044170240,      # NOT THE BEST KNOWN.  This is MW_128 (Elkies), the best known
                            # LATTICE kissing number at 128.  Edel-Rains-Sloane beat it by
                            # more than 40x in the same 1998 paper, at 8863556495104, which
                            # is in ERS_PUBLISHED.  Anyone extending floor_for past 96 who
                            # reaches for this entry will be low by a factor of 40.
}

# The two claims above that can be checked against each other, checked.
assert OTHER[80] > ERS_PUBLISHED[80], (OTHER[80], ERS_PUBLISHED[80])
assert OTHER[128] < ERS_PUBLISHED[128], (OTHER[128], ERS_PUBLISHED[128])
assert OTHER[64] == ERS_PUBLISHED[64], (OTHER[64], ERS_PUBLISHED[64])


# ---------------------------------------------------------------------------------------
# The UPPER-bound column of the same table, parsed from the cached copy at
# verifications/closed/dim32-44-ers-audit/data/cohn_kissing.html.  Dimensions 1, 2, 3, 4, 8
# and 24 are exact, so the upper bound equals the lower one there.  Nothing in this
# repository is DERIVED from these numbers; they exist so that audit.py can check that no
# claim exceeds a PROVEN upper bound, which is the one way a lower bound can be flatly
# impossible rather than merely unsupported.  Above dimension 48 the table has only n = 72,
# and monotonicity carries it down: tau(n) <= tau(72) <= 2545617287927 for every n <= 72.
COHN_UPPER = {
    1: 2, 2: 6, 3: 12, 4: 24, 5: 44, 6: 77,
    7: 134, 8: 240, 9: 363, 10: 553, 11: 868, 12: 1355,
    13: 2064, 14: 3174, 15: 4853, 16: 7320, 17: 10978, 18: 16406,
    19: 24417, 20: 36195, 21: 53524, 22: 80810, 23: 122351, 24: 196560,
    25: 265006, 26: 367775, 27: 522212, 28: 752292, 29: 1075991, 30: 1537707,
    31: 2213487, 32: 3162316, 33: 4494570, 34: 6422593, 35: 9162403, 36: 13017098,
    37: 18498316, 38: 26496684, 39: 37826766, 40: 53589200, 41: 76287040, 42: 108404055,
    43: 153813582, 44: 220788272, 45: 316735249, 46: 441900184, 47: 621658419, 48: 867897072,
    72: 2545617287927,
}


def upper_for(d):
    """(value, explanation) for the best proven upper bound applying in dimension d, or None"""
    if d in COHN_UPPER:
        return COHN_UPPER[d], "Cohn's table, dimension %d" % d
    above = [m for m in COHN_UPPER if m > d]
    if above:
        m = min(above)
        return COHN_UPPER[m], ("tau(%d) <= tau(%d) <= %d, Cohn's table at %d"
                               % (d, m, COHN_UPPER[m], m))
    return None


def floor_for(d):
    """(value, explanation) for the bound a claim in dimension d must beat.

    Published wherever a published entry exists: Cohn at 1-48 and 72, and the direct sums
    built on those.  Nothing published has an entry at 96: Cohn stops at 48 and 72, and the
    Nebe-Sloane table it replaced listed 1-40, 42, 44, 48, 64, 72, 80 and 128.  ERS 1998 is
    titled "in dimensions 32 to 128", but that is the RANGE OF A CONSTRUCTION, not a table of
    values: section 3 says "We illustrate the construction by giving new records in dimensions
    32, 36, 40, 44, 64, 80 and 128", and those seven are its only worked examples (checked
    against arXiv:math/0207291 on 2026-08-26).  So ERS tabulated neither 96 nor 62 nor 63.
    Their n = 64 and n = 80 totals, 331737984 and 1368532064, are the two reproduced in
    closed/dim32-44-ers-audit/ and both match the source to the digit.  Substituting current
    code tables into their formula is explicitly sanctioned by their Remark (1): any available
    lower bounds on A(n,d) and A(n,d,w) still give a lower bound on tau_n.

    Dimension 96 is therefore the same situation as 62 and 63 -- the ERS construction
    evaluated at a length with no record of its having been evaluated -- and is treated the
    same way: the floor is the direct sum Gamma_72 (+) Leech, and 12886999232 is a CLAIM.
    Until 2026-08-26 that value was the floor here instead, which removed 96 from the range by
    definition: "the range ends at 95" and "the floor at 96 is 12886999232" were one decision
    stated twice, not two facts, and the improvement it would otherwise have measured could
    not appear because the thing it improved on WAS it.  Beware that shape generally -- a
    baseline set to your own new result silently converts a claim into a non-result.
    """
    if d in COHN:
        return COHN[d], "Cohn's table: %s" % COHN_SOURCE.get(d, "see the table")
    if 49 <= d <= 63:
        return (COHN[48] + TAU[d - 48],
                "direct sum P_48 (+) best %d-dimensional lattice = 52416000 + tau(%d)"
                % (d - 48, d - 48))
    if 64 <= d <= 71:
        return (OTHER[64],
                "inherited by monotonicity from dimension 64 (Edel-Rains-Sloane 1998); "
                "no table has an entry for dimensions 65-71")
    if 73 <= d <= 95:
        # d = 80 is the one dimension in this range a table covers: Nebe-Sloane
        # record Gamma_72 perp E_8, and it is this same direct sum.  Say so, or the
        # generated provenance claims less than the repository can support.
        note = ("; also the Nebe-Sloane table entry Gamma_72 perp E_8 -- the only entry any "
                 "table has anywhere in 49-63, 65-71 or 73-96") if d in TABULATED_ABOVE_48 else ""
        return (COHN[72] + TAU[d - 72],
                "direct sum Gamma_72 (+) best %d-dimensional lattice = 6218175600 + tau(%d)%s"
                % (d - 72, d - 72, note))
    if d == 96:
        # COHN[24], not TAU[24]: TAU stops at 23, which is all dimensions 73-95 need.
        return (COHN[72] + COHN[24],
                "direct sum Gamma_72 (+) Leech = 6218175600 + 196560.  The ERS chain "
                "(96,24,6,1) reaches 12886999232 here, which is the CLAIM and not the floor "
                "-- see this function's docstring.  closed/dim96-ers-takeover/")
    raise ValueError(
        "no floor recorded for dimension %d.  Above 96 the Edel-Rains-Sloane count is the "
        "thing to evaluate and it has not been evaluated here; monotonicity gives at least "
        "dimension 96's %d, and nothing in this repository claims anything there."
        % (d, OTHER[96]))


# The one dimension a table covers inside 49-63, 65-71, 73-96.  If this ever grows, the
# sentences in sections 1.2, 1.5 and 8.2 of the paper that name 80 by itself are wrong.
_RANGES = list(range(49, 64)) + list(range(65, 72)) + list(range(73, 97))
TABULATED_ABOVE_48 = tuple(d for d in _RANGES if d in NEBE_SLOANE or d in COHN)
assert TABULATED_ABOVE_48 == (80,), TABULATED_ABOVE_48
assert floor_for(80)[0] == OTHER[80] == COHN[72] + 240, floor_for(80)

FLOOR = {d: floor_for(d)[0] for d in list(range(1, 97))}

REFERENCES = """
[Bia26]  F. Bianchi, Y. Kwon, A. Pappu, J. Zou, Harnessing the collective intelligence of
         AI agents in the wild for new discoveries, arXiv:2606.10402 (2026).
[BC23]   P. Boyvalenkov, D. Cherkashin, Results Math. 80 (2025), no. 1, Paper No. 3,
         doi:10.1007/s00025-024-02322-0; preprint arXiv:2312.05121 (Dec 2023), whose
         numbering the equation number below is.  That is the version read here, and
         the two need not number alike.  Equation (3) is the
         one-point distribution of P_48's minimal shell, and the same paragraph observes
         that its A_0 = 23766960 vectors form a 47-dimensional kissing configuration.
[Bro]    A. E. Brouwer, Bounds for binary constant weight codes,
         https://aeb.win.tue.nl/codes/Andw.html
[CL24]   H. Cohn, A. Li, Improved kissing numbers in seventeen through twenty-one
         dimensions, arXiv:2411.04916 (2024).
[Cohn]   H. Cohn, Table of kissing number bounds, https://cohn.mit.edu/kissing-numbers,
         https://hdl.handle.net/1721.1/153312
[CS99]   J. H. Conway, N. J. A. Sloane, Sphere Packings, Lattices and Groups, 3rd ed.,
         Springer 1999.
[Ech26]  W. Echols, New lower bounds for constant-weight codes via seeded bit-swap tabu
         search, arXiv:2608.13906 (2026).
[EdRS98] Y. Edel, E. M. Rains, N. J. A. Sloane, On kissing numbers in dimensions 32 to
         128, Electron. J. Combin. 5 (1998) R22, doi:10.37236/1360, arXiv:math/0207291.
[Ho26]   B. S. Ho, A new lower bound for the kissing number in 19 dimensions,
         arXiv:2603.10425 (2026).
[KKW16]  A. Kallal, D. Kan, E. Wang, SIAM J. Discrete Math. 31 (2017) 1895-1908,
         doi:10.1137/16M1095810; preprint arXiv:1608.07270 (2016), the year the
         chronologies elsewhere key this label to.  488 vectors in the class
         problem, the predecessor of PackingStar's 496.
[Lee67]  J. Leech, Notes on sphere packings, Canadian J. Math. 19 (1967) 251-267.
[LS71]   J. Leech, N. J. A. Sloane, Sphere packings and error-correcting codes, Canadian
         J. Math. 23 (1971) 718-745.
[Ma25]   C. Ma, T. Tao Z., P. Li, M. Liu, H. Chen, Z. Mao, Y. Cheng, Y. Qi, Y. Yang,
         Finding kissing numbers with game-theoretic reinforcement learning,
         arXiv:2511.13391 (2025).  "PackingStar"; holds dimensions 25-31.
[Neb12]  G. Nebe, An even unimodular 72-dimensional lattice of minimum 8, J. Reine Angew.
         Math. 673 (2012) 237-247, arXiv:1008.2862.
[NebCat] G. Nebe, N. J. A. Sloane, A Catalogue of Lattices,
         https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/  (entries P48p, P48q,
         Gamma72: Gram matrices and automorphism generators).
[Oze16]  M. Ozeki, A numerical study of Siegel theta series of various degrees for the
         48-dimensional even unimodular extremal lattices, Tsukuba J. Math. 40 (2016)
         139-186.
[SW26]   X. Sun, C. Wang, Sphere packings and kissing numbers in dimensions 39, 43 and 45
         from the antipode construction, arXiv:2607.20359 (2026).
[Tak26]  R. Takhanov, Z. Assylbekov, S. Yun, Structure of kissing arrangements in R^12 and
         a place for the 841st sphere, arXiv:2606.18984 (2026).
[Ven84]  B. B. Venkov, Even unimodular extremal lattices, Trudy Mat. Inst. Steklov 165
         (1984) 43-48.  The minimal shell of an extremal even unimodular lattice of
         dimension 0 mod 24 is a spherical 11-design.
[ZE99]   V. A. Zinoviev, T. Ericson, New lower bounds for contact numbers in small
         dimensions, Problems Inform. Transmission 35 (1999) 287-294.
[ZL84]   V. A. Zinoviev, S. N. Litsyn, On shortening of codes, Problems Inform.
         Transmission 20:1 (1984) 1-7.
"""


if __name__ == '__main__':
    print(__doc__)
    print("=" * 86)
    print("  dim   previously published   provenance")
    print("=" * 86)
    for d in range(1, 97):
        try:
            v, why = floor_for(d)
        except ValueError:
            continue
        print("  %3d   %-20d   %s" % (d, v, why))
    print()
    print("=" * 86)
    print("REFERENCES")
    print("=" * 86)
    print(REFERENCES)
