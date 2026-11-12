"""Exact obstruction for several cheap additional-height 37D candidate shapes.

This retains the complete current base and existing +/-2 latitude layer.
It does not cover arbitrary real bodies, altered base ownership, or general
short-shell patterns outside the explicitly screened shapes.
"""
from pathlib import Path
import json
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parents[2]/'verifications/improved/dim36-37-coxeter-todd'
sys.path.insert(0,str(SOURCE/'phase'))
from exact_lattice import enumerate_short


def main():
    c=json.loads((SOURCE/'completion/certificate.json').read_text())
    G=np.array(c['gram_12'],np.int64);M=np.array(c['sublattice_matrix'],np.int64)
    phi=np.eye(12,dtype=np.int64)-M
    R0=np.array(c['root_phase_zero_vectors'],np.int64)
    R=np.stack([R0,R0@phi,R0@phi@phi],axis=1)
    roots=R.reshape(-1,12)
    short=list(enumerate_short(G.tolist()))
    V=np.array([v for v,n in short],np.int64);N=np.array([n for v,n in short],np.int64)
    first,second=np.array(c['label_first']),np.array(c['label_second']);powers=3**np.arange(6)
    labels=((V[:,first]+V[:,second])%3)@powers
    weights=np.full(729,999,np.int64);weights[0]=0;np.minimum.at(weights,labels,N)
    six=V[N==6];sixlabels=labels[N==6];H=np.flatnonzero(weights==6)
    raw=(roots@G)@six.T
    maxima=np.array([raw[:,sixlabels==h].max(axis=1)for h in H]).T
    minimum=int(maxima.min());assert minimum>0
    base=np.array(c['base_outer_cosets'],np.int64)
    pair_coverage={}
    for k in range(3):
        rows=base[(base[:,k]==0)&np.all(np.sort(weights[base],axis=1)==[0,6,6],axis=1)]
        i,j=[q for q in range(3)if q!=k]
        assert len(rows)==224
        assert set(rows[:,i])==set(H)==set(rows[:,j])
        pair_coverage[str(k)]=len(rows)
    # Given body (x4,y6,0), choose the unchanged (6,6,0) owner whose second
    # coset contains y6. It contains the exact y6 and an x6 with x4 dot x6>0.
    # Its product with the proposed body is therefore strictly greater than six.
    lift=np.array(c['lift_phase_orbits'],np.int64)
    single={}
    for j in range(3):
        used=sorted(set(int(a)for a in lift[:,j]if a!=252))
        occupied=R[used].reshape(-1,12)
        fourmax=((roots@G)@occupied.T).max(axis=1)
        sixmax=((six@G)@occupied.T).max(axis=1)
        assert int(fourmax.min())>=1 and int(sixmax.min())>=2
        single[str(j)]={'used_root_labels':len(used),
                        'minimum_largest_product_for_norm4_body':int(fourmax.min()),
                        'minimum_largest_product_for_norm6_body':int(sixmax.min())}
        ownerlabels=set(int(row[j])for row in base if sorted(weights[row])==[0,4,8]and weights[row[j]]==8)
        assert ownerlabels==set(np.flatnonzero(weights==8))
    report={'status':'EXACT REJECTION OF SPECIFIED SHAPES',
            'norm4_plus_norm6_body_candidates':6*len(roots)*len(six),
            'unchanged_norm6_pair_owners_by_zero_block':pair_coverage,
            'minimum_blocking_base_product':6+minimum,
            'single_norm8_body_candidates':3*int(np.sum(N==8)),
            'single_norm8_blocking_base_product':8,
            'single_block_norm4_and_norm6':single,
            'head_squares':{'body_norm4':8,'body_norm6':6,'body_norm8':4,'body_norm10':2},
            'reason':'Norm10 mixed bodies and norm8 single bodies already conflict with the unchanged equatorial base. Norm4/6 single bodies conflict with one of the retained +/-2 latitude points.',
            'scope':'No changes to base ownership or existing latitude bodies; only the explicitly enumerated integer lattice shapes.'}
    (HERE/'other37_latitudes_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
