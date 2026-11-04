"""Enumerate the two saturated 2+2 bridge-coupling types; see BRIDGES.md."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import hashlib,itertools,json
from pathlib import Path
import numpy as np
from alternative_bridges import bridge_supports
ROOT=Path(__file__).resolve().parent


def main():
    _,_,matchings,factor_sets=bridge_supports();edges=list(itertools.combinations(range(6),2))
    edge_indices=[[edges.index(tuple(e)) for e in m] for m in matchings]
    factors=[sorted(f) for f in factor_sets];standard=set();alternative=set()
    for left,right in itertools.product(factors,repeat=2):
        for order in itertools.permutations(right):
            matrix=np.zeros((15,15),dtype=np.uint8)
            for ml,mr in zip(left,order):matrix[np.ix_(edge_indices[ml],edge_indices[mr])]=1
            assert np.all(matrix.sum(axis=0)==3) and np.all(matrix.sum(axis=1)==3)
            standard.add(np.packbits(matrix.ravel()).tobytes())
    for order in itertools.permutations(range(6)):
        matrix=np.zeros((15,15),dtype=np.uint8)
        for i,(u,v) in enumerate(edges):
            common=factor_sets[order[u]]&factor_sets[order[v]];assert len(common)==1
            matrix[i,edge_indices[next(iter(common))]]=1
        alternative.add(np.packbits(matrix.ravel()).tobytes())
    assert len(standard)==4320 and len(alternative)==720 and not standard&alternative
    # Verify that all 5! permutations of one factorization's colors are induced
    # by vertex relabelings, as used in the two-type classification.
    matching_lookup={tuple(sorted(tuple(e) for e in m)):i for i,m in enumerate(matchings)}
    color_actions=set();factor_orbit=set();fixed=factors[0]
    for p in itertools.permutations(range(6)):
        image=[]
        for m in fixed:
            transformed=tuple(sorted(tuple(sorted((p[a],p[b]))) for a,b in matchings[m]))
            image.append(matching_lookup[transformed])
        factor_orbit.add(tuple(sorted(image)))
        if set(image)==set(fixed):color_actions.add(tuple(fixed.index(v) for v in image))
    assert len(color_actions)==120 and len(factor_orbit)==6
    report=dict(exact_enumeration=True,standard_labeled_couplings=len(standard),
                alternative_labeled_couplings=len(alternative),total_labeled_couplings=5040,
                standard_sha256=hashlib.sha256(b''.join(sorted(standard))).hexdigest(),
                alternative_sha256=hashlib.sha256(b''.join(sorted(alternative))).hexdigest(),
                all_120_color_permutations_realized=True,all_six_factorizations_isomorphic=True,
                scope='45 distinct 2+2 supports, all sixteen signs on each support; see proof in BRIDGES.md')
    (ROOT/'couplings.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
