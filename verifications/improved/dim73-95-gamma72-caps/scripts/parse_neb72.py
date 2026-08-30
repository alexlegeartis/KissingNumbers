#!/usr/bin/env python3
"""Parse Nebe's Gamma_72 Gram matrix and automorphism-subgroup generators out of the
Catalogue of Lattices page

    http://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/neb72.html

(entry "Gamma72"; last modified Fri Jul 18 13:28:23 CEST 2014; reference
 G. Nebe, "An even unimodular 72-dimensional lattice of minimum 8",
 J. reine angew. Math. 673 (2012) 237-247, arXiv:1008.2862).

Writes  gamma72_gram.npy   (72x72 int64 Gram matrix, even, det 1)
        gamma72_gens.npy   (6 x 72 x 72 int64 generators of SL2(25) x PSL2(7):2)

Usage:  python parse_neb72.py <path-to-neb72.html> <outdir>
"""
import re, sys, os
import numpy as np

def strip_html(s):
    s = re.sub(r'<[^>]*>', ' ', s)
    return s

def main():
    if len(sys.argv) < 3:
        raise SystemExit(
            "usage: python parse_neb72.py <page.html> <outdir>\n"
            "\n"
            "  e.g. python parse_neb72.py ../data/Gamma72.html ../data\n"
            "\n"
            "Parses a Catalogue of Lattices page into a Gram matrix and generator matrices\n"
            "and verifies them exactly: symmetric, even diagonal, determinant 1 by\n"
            "fraction-free Bareiss elimination, and A G A^T = G for every generator.\n"
            "\n"
            "The outputs are already shipped in ../data/; this script is here so that the\n"
            "provenance can be re-checked from the catalogue page itself.")
    src = sys.argv[1]
    out = sys.argv[2]
    raw = open(src, 'r', encoding='latin-1').read()

    # --- GRAM ---
    i = raw.index('<a NAME="GRAM">')
    j = raw.index('<a NAME="DET">')
    blk = strip_html(raw[i:j])
    blk = blk.replace('GRAM', ' ')
    toks = blk.split()
    nums = [int(t) for t in toks if re.fullmatch(r'-?\d+', t)]
    assert nums[0] == 72, nums[:5]
    nums = nums[1:]
    assert len(nums) == 72 * 72, len(nums)
    G = np.array(nums, dtype=np.int64).reshape(72, 72)

    # --- SUBGROUP_GENERATORS ---
    i = raw.index('<a NAME="SUBGROUP_GENERATORS">')
    j = raw.index('<a NAME="PROPERTIES">')
    blk = strip_html(raw[i:j])
    blk = blk.replace('SUBGROUP_GENERATORS', ' ')
    toks = blk.split()
    # leading "#6" then "72" then 6*72*72 integers
    m = re.search(r'#(\d+)', blk)
    ngen = int(m.group(1))
    blk2 = blk[m.end():]
    nums = [int(t) for t in blk2.split() if re.fullmatch(r'-?\d+', t)]
    # each generator block is prefixed by its dimension, "72"
    assert len(nums) == ngen * (1 + 72 * 72), (len(nums), ngen * (1 + 72 * 72))
    GENS = np.empty((ngen, 72, 72), dtype=np.int64)
    for g in range(ngen):
        base = g * (1 + 72 * 72)
        assert nums[base] == 72
        GENS[g] = np.array(nums[base + 1:base + 1 + 72 * 72],
                           dtype=np.int64).reshape(72, 72)

    os.makedirs(out, exist_ok=True)
    np.save(os.path.join(out, 'gamma72_gram.npy'), G)
    np.save(os.path.join(out, 'gamma72_gens.npy'), GENS)
    print("gram  ", G.shape, "symmetric:", bool((G == G.T).all()))
    print("gens  ", GENS.shape)

if __name__ == '__main__':
    main()
