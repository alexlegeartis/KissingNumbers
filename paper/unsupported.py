"""Every distinct integer >= 1000 in kissing46.tex, checked to appear somewhere else in
the repository.  A figure in the paper with no home in the verification packages is
unsupported, and that is what this catches.  Two kinds of false positive are expected
and harmless: numbers the paper derives on the page (the Table 1 gains), and digits
glued across a comma-separated list such as "n = 1,2,3,4,8,24".

    python unsupported.py      # a second, no third-party packages
"""
import os, re, io

HERE = os.path.dirname(os.path.abspath(__file__))


def _find_repo():
    """The verification repository, whether or not this file sits inside it."""
    here = os.path.dirname(HERE)
    for cand in (here, os.path.join(here, 'kissing_verifications')):
        if os.path.exists(os.path.join(cand, 'RESULTS.md')):
            return cand
    raise SystemExit('cannot find the verification repository: expected RESULTS.md in\n'
                     '  %s\nor\n  %s' % (here, os.path.join(here, 'kissing_verifications')))


ROOT = _find_repo()
TEX = io.open(os.path.join(HERE, 'kissing46.tex'), encoding='utf-8').read()

SEP = re.compile(u'[\\s,\u00a0\u202f\\\\]+')


def digits(s):
    return SEP.sub('', s)


# collect every number in the paper, allowing LaTeX \, separators
nums = set()
for m in re.finditer(r'\d[\d\\, ]*\d|\d', TEX):
    t = digits(m.group(0))
    if t.isdigit() and int(t) >= 1000:
        nums.add(int(t))

corpus = {}
for root, ds, fs in os.walk(ROOT):
    ds[:] = [d for d in ds if d not in ('__pycache__', '.git')]
    for f in fs:
        if not f.endswith(('.py', '.md', '.tex', '.txt', '.cff', '.bib', '.json')):
            continue
        p = os.path.join(root, f)
        if p.endswith('kissing46.tex'):
            continue
        try:
            corpus[os.path.relpath(p, ROOT)] = digits(
                io.open(p, encoding='utf-8', errors='replace').read())
        except Exception:
            pass

missing = []
for v in sorted(nums):
    s = str(v)
    where = [k for k, txt in corpus.items() if s in txt]
    if not where:
        missing.append(v)

# Figures the paper derives on the page, and so does not need a home elsewhere: the
# Table 1 gains and the dimension 81 - dimension 80 difference, all of which are
# differences of two columns of RESULTS.md, and the value a 249-line class would give
# in dimension 25, which is a hypothetical, not a claim.
def _results_values():
    """{dimension: claimed value} from RESULTS.md, the file every claim is written from"""
    txt = io.open(os.path.join(ROOT, 'RESULTS.md'), encoding='utf-8').read()
    out = {}
    for m in re.finditer(r'^\| (\d+) \| ([0-9\u202f ]+) \| \*\*([0-9\u202f ]+)\*\*',
                         txt, re.M):
        out[int(m.group(1))] = (int(re.sub(r'[^0-9]', '', m.group(3))),
                                int(re.sub(r'[^0-9]', '', m.group(2))))
    return out


_PAIRS = _results_values()
_VAL = {d: v[0] for d, v in _PAIRS.items()}
# Table 1's gain column for "layers over Gamma_72" is the largest gain in 73-95, and it
# sat in DERIVED as a literal until the odd-k partitions moved it on 2026-08-30 -- the
# same fate the comment below records for the dimension-81 margin.  Computed now.
_G7395 = max((v - q for d, (v, q) in _PAIRS.items() if 73 <= d <= 95), default=None)
# The dimension 81 - dimension 80 margin is the difference of two rows of RESULTS.md, so
# it is computed rather than listed.  The
# margin used to sit here as a literal, with a comment recording that it had already moved
# once (184 904 -> 185 040); it moved twice more, on 2026-08-29, which is what a literal in
# a checker always ends up doing.
_MARGIN = _VAL[81] - _VAL[80] if 81 in _VAL and 80 in _VAL else None

DERIVED = {184980, 197060, 14947188, 2271920766,
           # what the icosahedral direction set would give in dimension 27, computed
           # from (4) in section 5.5 to show the cuboctahedron is forced
           199548,
           # the two Table 1 group maxima created when dimensions 62 and 63 moved out of
           # the P_48 group: 59898694 - 52417154 at dimension 61, and
           # 134217728 - 52418564 at dimension 63
           7481540, 86001280,
           # the Table 1 gain in dimension 96, 12886999232 - 6218372160
           6668627072,
           # the three levels below weight 96 of the dimension-96 chain, named in 6.3 as
           # the fixed part when the level-zero code is varied: 4^9*16384 + 16^4*32 + 96*2.
           # formulas.py asserts it.
           4297064640} | {v for v in (_MARGIN, _G7395) if v is not None}
# Artefacts of gluing digits across a separator, not numbers at all.
GLUED = {1234824, 184278}

print('distinct integers >= 1000 in the paper: %d' % len(nums))
print('of those, NOT found anywhere else in the repository: %d' % len(missing))
for v in missing:
    ctx = []
    for m in re.finditer(re.escape(str(v)[:4]), TEX):
        seg = TEX[max(0, m.start() - 90):m.start() + 60].replace('\n', ' ')
        if digits(TEX[m.start():m.start() + 30]).startswith(str(v)):
            ctx.append(re.sub(r'\s+', ' ', seg))
            break
    print('  %-14d %s' % (v, ctx[0] if ctx else ''))


unexplained = [v for v in missing if v not in DERIVED and v not in GLUED]
print()
print('=' * 78)
if unexplained:
    print('*** UNSUPPORTED FIGURES IN THE PAPER: %s ***'
          % ', '.join(str(v) for v in unexplained))
else:
    print('ALL CHECKS PASSED -- every figure in the paper is either derived on the page')
    print('or appears in a verification package.')
print('=' * 78)
raise SystemExit(1 if unexplained else 0)
