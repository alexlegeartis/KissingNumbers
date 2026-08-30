#!/usr/bin/env python3
"""Check every file of this package against data/SHA256SUMS.txt.

    python scripts/verify_manifest.py

This used to be documented as a manual step (`sha256sum -c data/SHA256SUMS.txt`), and
because nothing ran it, an edit to paper/kissing27.tex and a rebuild of the pdf once
invalidated the manifest without anything noticing. It is a test now.

Note for anyone editing this package: the manifest covers text files, so a line-ending
conversion invalidates it. The repository's .gitattributes disables end-of-line
conversion for exactly that reason.
"""
import hashlib
import io
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
MAN = os.path.join(PKG, 'data', 'SHA256SUMS.txt')


def main():
    if not os.path.exists(MAN):
        print("FAIL   no manifest at %s" % MAN)
        return 1

    entries = [l for l in io.open(MAN, encoding='utf-8').read().split('\n') if l.strip()]
    bad, missing = [], []
    for line in entries:
        want, rel = line.split('  ', 1)
        p = os.path.join(PKG, rel.replace('/', os.sep))
        if not os.path.exists(p):
            missing.append(rel)
            continue
        got = hashlib.sha256(open(p, 'rb').read()).hexdigest()
        if got != want:
            bad.append((rel, want, got))

    # anything tracked but absent from the manifest is a gap in the record
    tracked = {l.split('  ', 1)[1] for l in entries}
    extra = []
    for root, ds, fs in os.walk(PKG):
        ds[:] = [d for d in ds if d not in ('__pycache__', '.git')]
        for f in fs:
            rel = os.path.relpath(os.path.join(root, f), PKG).replace(os.sep, '/')
            if rel.startswith('data/SHA256SUMS') or rel.endswith('.pyc'):
                continue
            if rel not in tracked:
                extra.append(rel)

    print(__doc__.split('\n\n')[0])
    print("=" * 78)
    print("  manifest entries : %d" % len(entries))
    print("  verified         : %d" % (len(entries) - len(bad) - len(missing)))
    for rel in missing:
        print("  MISSING          : %s" % rel)
    for rel, want, got in bad:
        print("  MISMATCH         : %s" % rel)
        print("      expected %s" % want)
        print("      found    %s" % got)
    for rel in sorted(extra):
        print("  NOT IN MANIFEST  : %s" % rel)
    print("=" * 78)
    ok = not (bad or missing or extra)
    print(("ALL %d CHECKSUMS VERIFY" % len(entries)) if ok
          else "*** MANIFEST DOES NOT VERIFY ***")
    print("=" * 78)
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
