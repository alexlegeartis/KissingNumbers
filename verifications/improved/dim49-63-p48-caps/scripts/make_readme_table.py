#!/usr/bin/env python3
"""Rewrite the result table of README.md from the driver's own output.

    python final.py && python scripts/make_readme_table.py

The table had drifted badly: it was still on the superseded lambda = tau_lat(k)/2, so seven of
its fifteen rows were low and dimension 63 read 66 265 628 against the repository's 67 334 086.
Generating it removes the possibility.
"""
import io, json, os, re

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
DELIV = os.path.join(PKG, '..', '..', '..')
rows = json.load(open(os.path.join(PKG, 'rows49-63.json')))
readme = os.path.join(PKG, 'README.md')
text = io.open(readme, encoding='utf-8', newline='').read()

# Join with the README's own line ending: it is read with newline='' so the surrounding
# lines keep theirs, and a bare newline here would leave a CRLF file mixed.
EOL = '\r\n' if '\r\n' in text else '\n'
NB = u' '


def sp(n):
    return '{:,}'.format(int(n)).replace(',', NB)


# Dimensions 62 and 63 are computed but no longer claimed from this construction: a single
# Edel-Rains-Sloane sign code beats even its CEILING there (scripts/ceiling.py), so they moved
# to ../dim62-63-ers-chain/.  Dimension 61 is now the last row this package claims, so it is
# the one emphasised, and the two superseded rows say so in place of a bold value.
SUPERSEDED = {}
for _line in io.open(os.path.join(DELIV, 'RESULTS.md'), encoding='utf-8'):
    _m = re.match(r'\| (62|63) \| [^|]* \| \*\*([\d\s  ]+)\*\* \|', _line)
    if _m:
        SUPERSEDED[int(_m.group(1))] = int(re.sub(r'\D', '', _m.group(2)))
assert set(SUPERSEDED) == {62, 63}, 'RESULTS.md has no claim in dimension 62 or 63'
out = []
for r in rows:
    d = r['dim']
    b = (lambda x: '**%s**' % x) if d == 61 else (lambda x: x)
    if d in SUPERSEDED:
        val = '%s &nbsp;*(superseded by %s)*' % (sp(r['value']), sp(SUPERSEDED[d]))
    else:
        val = '**%s**' % sp(r['value'])
    out.append('| %s | %s | %s | %s | %s | %s | %s |'
               % (b(str(d)), b(str(r['k'])), b(sp(r['lam'])), b(sp(r['S'])),
                  b('+' + sp(r['gain'])), val, sp(r['previous'])))
new = EOL.join(out)

lines = text.split(EOL)
first = next(i for i, l in enumerate(lines) if l.startswith('| 49 |'))
last = max(i for i, l in enumerate(lines)
           if l.startswith('| 63 |') or l.startswith('| **63** |'))
text = EOL.join(lines[:first] + new.split(EOL) + lines[last + 1:])
io.open(readme, 'w', encoding='utf-8', newline='').write(text)
print("README.md rows 49-63 regenerated from rows49-63.json (%d rows)" % len(rows))
print(new.encode('ascii', 'replace').decode())
