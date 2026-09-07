# -*- coding: utf-8 -*-
"""Census of the prose faults listed in writing_guidelines.md.

    python style_scan.py [kissing46.tex] [-v]

Adapted from EDMGrokking/icomp_v2/style_scan.py. This paper has no appendix, so the
whole body is scanned as one segment; the theorem environments and the longtable are
stripped first, since displayed mathematics is not prose.

Run before and after every editing pass: passes reliably reintroduce what they fix.
Read the -v output; several checks are advisory and every hit needs a decision.
"""
import io
import os
import re
import sys

PATH = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith('-') \
    else os.path.join(os.path.dirname(os.path.abspath(__file__)), 'kissing46.tex')

SUBJ = (r"(?:it|they|we|this|that|these|those|there|the|its|their|each|every|both|neither|either|"
        r"no|one|two|three|four|five|six|seven|eight|nine|ten|all|most|some|nothing|only)")
VERB = (r"(?:is|are|was|were|has|have|had|does|do|did|can|could|may|might|must|will|would|"
        r"\w+(?:s|es|ed))")
NUM = r"(?:two|three|four|five|six|seven|eight|nine|ten)"

BARE_NUM = (r"\bthe " + NUM + r"\b(?=[,.;:]|\s+(?:is|are|was|were|do|does|did|fail|fails|failed|"
            r"show|shows|give|gives|report|reports|and|but|which|that|however|above|below)\b)")

COUNT_NOUN = (r"(?:statistics|things|observations|qualifications|reasons|properties|departures|"
              r"checks|results|consequences|conditions|features|points|notes|arguments|"
              r"explanations|comparisons|criteria|caveats|failures|ways|kinds|respects|"
              r"mechanisms|constraints|remarks)")
ENUM = (r"(?:First|Second|Third|The first|The second|One of|\\emph\{|\\item|\\textbf\{|"
        r"\\paragraph\{|Its |In the first|: )")
COUNT_HEAD = (r"\b(?:Two|Three|Four|Five|Six|Seven|Eight|Nine)\s+(?:further\s+|other\s+|more\s+)?"
              + COUNT_NOUN + r"\b(?![^.]*\bof\b)(?![\s\S]{0,240}?" + ENUM + r")")

OBLIQUE = (r"\bthe (?:one|only|single) (?:[a-z]+ ){0,2}"
           r"(?:candidate|alternative|option|choice|example|instance|case) (?:the|that|which)\b")

# Metaphor and informality where a plain verb exists. The list from the source guidelines,
# plus the words this manuscript actually reached for.
BANNED = (r"\b(?:sits?|sitting|buys?|trades?|exchanges?|cheap|cheaply|knobs?|arms?\b|pipelines?|"
          r"bites?\b|kills?|killed|killing|pays?\b|paid\b|earns?|dies?\b|died|die\b|"
          r"tempting|pretty|prettiest|dangerous|silently wrong|the whole game|"
          r"lottery|room\b|doors?\b|rescues?|manufactures?|flatters?|fires\b|the ladder|"
          r"absorbed|worth having|got wrong|goes? wrong|went wrong|nobody|anyone|"
          r"stale|leverage|good enough|cannot see|does not see|knows nothing)\b")

WHAT = r"\bwhat\b"
IS_WHAT = r"\b(?:is|are|was|were)\s+what\b"
PERIPHRASIS = (r"\b(?:is|are)\s+the\s+(?:quantity|one|thing|reason|value|case|version|point)\s+"
               r"(?:that|which)\b")

BARE_PAIR = (r"\bthe (?:two|three|four) (?:settings?|diagnostics|groups?|conditions?|statistics|"
             r"quantities|extremes|removals|failures|candidates|halves|arms|regimes|"
             r"rules|norms|spaces|sweeps|controls|findings|qualifications|departures|checks|"
             r"exceptions|measurements|mechanisms|constructions|families|ranges|routes)\b")

CONTRAST = r"\brather than\b|\bnot a\b|\bnot as\b|\bnone of\b|\band not\b"

BARE_DET = (r"(?:Both|Neither|Either|All (?:three|four|five|six)|The (?:first|last|other) "
            + NUM + r")")
FINITE = (r"(?:is|are|was|were|do|does|did|lie|lies|show|shows|report|reports|fail|fails|fix|"
          r"fixes|ask|asks|hold|holds|matter|matters|remain|remains|apply|applies|give|gives|"
          r"carr(?:y|ies)|separat\w+|explain\w*|behav\w+|appear\w*|say|says|come|comes)\b")
COUNT_SUBJ = r"(?:^|(?<=[.;:] ))" + BARE_DET + r"\s+" + FINITE

# Announcing what the next sentences will do.
ANNOUNCE = (r"\b(?:Two|Three|Four|Five) (?:features|things|consequences|remarks|points|"
            r"observations|facts)\b|\bare worth (?:recording|stating|noting)\b|"
            r"\bwe (?:record|note|remark|observe|emphasise|emphasize) (?:that|here)\b|"
            r"\bis worth (?:recording|stating|noting)\b|\bIt is worth\b")

# The author instructing the author, or narrating the manuscript's own history.
SELFTALK = (r"\b(?:an earlier version|earlier version of (?:our|this)|our own code|"
            r"we recommend|we would (?:most )?like|we say so|we flag|we claim neither|"
            r"we do not count|we have not found|we could not|we retain|had to be|"
            r"which is why|the symptom|we present them as such)\b")

CHECKS = [
    ("pseudo-cleft",           r"\bWhat [a-z][^.]{0,60}? is\b", True),
    ("'where' for 'whereas'",  r"[a-z], where [a-z]", True),
    ("initial And/But",        r"\. (?:And|But) ", True),
    ("inanimate 's",           r"\b[a-z][a-z-]+'s\b", True),
    ("banned word",            BANNED, True),
    ("', and ' clause-join",   r", and (" + SUBJ + r"\b[^,.;:]{0,60}?\b" + VERB + r"\b)", True),
    ("em dash",                r"---", False),
    ("'what' clause",          WHAT, True),
    ("'is what'",              IS_WHAT, True),
    ("'is the X that'",        PERIPHRASIS, True),
    ("count as subject",       COUNT_SUBJ, True),
    ("bare pair",              BARE_PAIR, True),
    ("bare numeral for a set", BARE_NUM, True),
    ("count, items not listed", COUNT_HEAD, True),
    ("oblique naming",         OBLIQUE, True),
    ("announcing the next",    ANNOUNCE, True),
    ("self-talk / history",    SELFTALK, True),
    ("contrast (advisory)",    CONTRAST, True),
    ("', which'",              r", which ", False),
    ("', so '",                r", so ", False),
    ("first person 'we'",      r"\bwe\b|\bour\b|\bus\b", False),
    ("third person author",    r"\bthe (?:present )?author(?:'s)?\b", True),
]


def flatten(s):
    for env in ('tabular', 'longtable', 'align', 'align\\*', 'equation', 'equation\\*',
                'array', 'itemize', 'center'):
        s = re.sub(r"\\begin\{%s\}.*?\\end\{%s\}" % (env, env), " ", s, flags=re.S)
    s = re.sub(r"\\\[.*?\\\]", " EQ ", s, flags=re.S)
    s = re.sub(r"\$[^$]*\$", " X ", s)
    s = re.sub(r"\\(?:cite|ref|eqref|label)\{[^}]*\}", " ", s)
    return re.sub(r"\s+", " ", s)


def sentences(flat):
    return [t for t in re.split(r"(?<=[.!?]) (?=[A-Z\\])", flat) if t.strip()]


def main():
    text = io.open(PATH, encoding='utf8').read()
    body = text[text.index(r"\begin{abstract}"):text.index(r"\bibliographystyle")]
    flat = flatten(body)
    words = len(flat.split())
    sents = sentences(flat)
    verbose = '-v' in sys.argv

    print("== %s  %d words, %d sentences, %.1f words/sentence"
          % (os.path.basename(PATH), words, len(sents), words / len(sents)))
    for label, pat, ctx in CHECKS:
        hits = list(re.finditer(pat, flat, re.I))
        print("   %-26s %4d   (%.1f / 1000 words)"
              % (label, len(hits), 1000.0 * len(hits) / words))
        if verbose and ctx:
            for m in hits[:14]:
                print("        ...%s" % flat[max(0, m.start() - 66):m.end() + 60])

    long_s = [s for s in sents if len(s.split()) > 38]
    print("   %-26s %4d" % ("sentences > 38 words", len(long_s)))
    if verbose:
        for s in long_s[:10]:
            print("        (%d) %s" % (len(s.split()), s[:190]))

    caps = []
    for m in re.finditer(re.escape("\\caption{"), text):
        i, d, j = m.end(), 1, m.end()
        while d:
            d += (text[j] == '{') - (text[j] == '}')
            j += 1
        lab = re.search(re.escape("\\label{") + r"([^}]+)\}", text[j:j + 160])
        caps.append((lab.group(1) if lab else '?', len(flatten(text[i:j - 1]).split())))
    if caps:
        print("== CAPTIONS  %d, mean %.0f words"
              % (len(caps), sum(n for _, n in caps) / len(caps)))
        for lab, n in sorted(caps, key=lambda c: -c[1]):
            print("   %4d  %s" % (n, lab))


main()
