"""Expose 756 independent fixed-base choices and retain compatible seed orbits.

This constructs search inputs, not an improved bound. Both complete exact
verifiers must check any final certificate after a phase/lift refill.
"""
import argparse
import hashlib
import json
from pathlib import Path
import random
import sys
import numpy as np

PACKAGE=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(PACKAGE/'phase'))
from exact_lattice import enumerate_short


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificate',type=Path,default=PACKAGE/'completion/certificate.json')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--flips',type=int,default=1)
    parser.add_argument('--seed',type=int,default=362810)
    args=parser.parse_args()
    output=args.output.resolve()
    if output==PACKAGE or PACKAGE in output.parents:
        parser.error('write proposals outside the verification package')
    raw=args.certificate.read_bytes();c=json.loads(raw)
    G=np.array(c['gram_12'],dtype=np.int64)
    powers=3**np.arange(6,dtype=np.int64)
    weights=np.full(729,100,dtype=np.int16);weights[0]=0
    for v,n in enumerate_short(G.tolist()):
        label=sum(((v[a]+v[b])%3)*p for a,b,p in zip(c['label_first'],c['label_second'],powers))
        weights[label]=min(weights[label],n)
    digits=np.arange(729)[:,None]//powers%3
    D=weights[((digits[:,None,:]-digits[None,:,:])%3)@powers]
    negative=((-digits)%3)@powers
    roots=np.flatnonzero(weights==4)
    base=np.array(c['base_outer_cosets'],dtype=np.int16)
    assert base.shape==(1429,3)
    alternatives=[];owners=[]
    for index,row in enumerate(base):
        norms=weights[row]
        if sorted(map(int,norms))!=[0,4,8]:
            continue
        i=int(np.flatnonzero(norms==4)[0]);j=int(np.flatnonzero(norms==8)[0])
        k=3-i-j
        alternative=[0,0,0]
        alternative[j]=int(row[j]);alternative[k]=int(negative[row[i]])
        owners.append(index);alternatives.append(alternative)
    A=np.array(alternatives,dtype=np.int16)
    assert A.shape==(756,3) and len(set(owners))==756
    assert len(set(map(tuple,A.tolist())))==756
    # The zero coset represents the L^3 shell, handled separately in the proof.
    original=sum(D[base[:,j,None],base[:,j]] for j in range(3))
    np.fill_diagonal(original,12)
    assert int(original.min())>=12
    against_old=sum(D[A[:,j,None],base[:,j]] for j in range(3))
    against_old[np.arange(756),owners]=12
    assert int(against_old.min())>=12
    against_new=sum(D[A[:,j,None],A[:,j]] for j in range(3))
    np.fill_diagonal(against_new,12)
    assert int(against_new.min())>=12
    # Thus every subset of these 756 replacements retains all base angles.
    if not 0<=args.flips<=756:
        parser.error('--flips must be between 0 and 756')
    selected=sorted(random.Random(args.seed).sample(range(756),args.flips))
    proposed=base.copy()
    for index in selected:
        proposed[owners[index]]=A[index]
    c['base_outer_cosets']=proposed.tolist()
    phase=[]
    for orbit in c['phase_orbits']:
        labels=roots[orbit[:3]]
        if np.all(sum(D[labels[j],proposed[:,j]] for j in range(3))>=12):
            phase.append(orbit)
    lifts=[]
    for orbit in c['lift_phase_orbits']:
        labels=[0 if a==252 else int(roots[a]) for a in orbit[:3]]
        if np.all(sum(D[labels[j],proposed[:,j]] for j in range(3))>=8):
            lifts.append(orbit)
    c['phase_orbits']=phase;c['lift_phase_orbits']=lifts
    c['points_36']=403704+3*len(phase)
    c['points_37']=c['points_36']+6*len(lifts)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(c,separators=(',',':'))+'\n')
    print(json.dumps(dict(base_alternatives=756,every_subset_preserves_base=True,
        base_points=403704,flipped_indices=selected,seed=args.seed,
        retained_phase_orbits=len(phase),retained_lift_orbits=len(lifts),
        points_36=c['points_36'],points_37=c['points_37'],
        source_sha256=hashlib.sha256(raw).hexdigest(),
        scope='Search input; no improvement or completed refill is claimed.')))


if __name__=='__main__':
    main()
