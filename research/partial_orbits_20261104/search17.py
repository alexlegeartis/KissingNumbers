"""Construct a 17D seed and search individual phases, with no fixed layers."""
import argparse
from itertools import combinations
import json
from pathlib import Path
import sys
from time import monotonic

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from discovery.core import Artifact,RunContext,digest
from discovery.certification import verify_artifact
from discovery.engines.partial_orbits import select_partial_orbits
from discovery.runner import source_fingerprints


def binary_span_basis(words):
    pivots={}
    for word in words:
        while word:
            p=word.bit_length()-1
            if p in pivots:word^=pivots[p]
            else:pivots[p]=word;break
    return list(pivots.values())


def sign_reps(positions,group):
    projection=[sum(((h>>j)&1)<<i for i,j in enumerate(positions)) for h in group]
    basis=binary_span_basis(projection)
    pivots={b.bit_length()-1 for b in basis}
    free=[i for i in range(len(positions)) if i not in pivots]
    return [sum(((w>>i)&1)<<j for i,j in enumerate(free)) for w in range(1<<len(free))]


def seed_and_group():
    octads=json.loads((ROOT/'research/mixed_tier17_20261104/octads18.json').read_text())['octads']
    masks=[sum(1<<i for i in o) for o in octads]
    group=sorted(set([0,65535]+masks))
    assert len(group)==32
    rm2=[w for w in range(65536) if all((w&o).bit_count()%2==0 for o in masks)]
    remaining=set(rm2);cosets=[]
    while remaining:
        w=min(remaining);cosets.append(w);remaining.difference_update(w^g for g in group)
    assert len(cosets)==64
    signs=np.array([[1-2*((g>>i)&1) for i in range(16)] for g in group],np.int64)
    plus=[0,12,22,27,39,41];minus=[1,13,23,26,38,40]
    flats={h:np.array([[1-2*(((cosets[c]^g)>>i)&1) for i in range(16)]
                     for c in indices for g in group],np.int64)
           for h,indices in [(4,plus),(-4,minus)]}
    short=[]
    for support in (w for w in rm2 if w.bit_count()==6):
        for w in (0,support&-support):
            v=np.array([(1-2*((w>>i)&1))*2 if support>>i&1 else 0 for i in range(16)])
            orbit=v*signs
            if np.max(orbit@flats[4].T)<=8 and np.max(orbit@flats[-4].T)<=8:
                short.extend(orbit.tolist())
    assert len(short)==512
    rows=[]
    for a,b in combinations(range(16),2):
        for x in (-12,12):
            for y in (-12,12):
                row=[0]*17;row[a]=x;row[b]=y;rows.append(row)
    for support in masks:
        positions=[i for i in range(16) if support>>i&1]
        for w in range(256):
            if w.bit_count()%2:
                row=[0]*17
                for j,i in enumerate(positions):row[i]=6*(1-2*((w>>j)&1))
                rows.append(row)
    for h in (-6,6):rows.extend([[3*x for x in v]+[h] for v in short])
    for h,values in flats.items():rows.extend([[4*int(x) for x in v]+[h] for v in values])
    rows.extend([[0]*16+[12],[0]*16+[-12]])
    assert len(rows)==5730 and len(set(map(tuple,rows)))==5730
    basis=binary_span_basis(group)
    generators=[[(g>>i)&1 for i in range(16)]+[0] for g in basis]
    return rows,group,generators,masks


def representatives(group,masks,families):
    rows=[]
    if 'six' in families:
        for positions in combinations(range(16),6):
            for w in sign_reps(positions,group):
                body=[0]*16
                for j,i in enumerate(positions):body[i]=6*(1-2*((w>>j)&1))
                rows.extend([body+[6],body+[-6]])
    if 'dual' in families:
        for q in range(16):
            for w in sign_reps(list(range(16)),group):
                body=[(9 if i==q else 3)*(1-2*((w>>i)&1)) for i in range(16)]
                rows.extend([body+[6],body+[-6]])
    if 'flat' in families:
        for w in sign_reps(list(range(16)),group):
            body=[4*(1-2*((w>>i)&1)) for i in range(16)]
            rows.extend([body+[4],body+[-4]])
    if 'newshort' in families:
        for q in range(16):
            for others in combinations([i for i in range(16) if i!=q],2):
                positions=[q,*others]
                for w in sign_reps(positions,group):
                    body=[0]*16
                    for j,i in enumerate(positions):body[i]=(12 if j==0 else 6)*(1-2*((w>>j)&1))
                    rows.extend([body+[6],body+[-6]])
    if 'equator' in families:
        for positions in combinations(range(16),8):
            for w in sign_reps(positions,group):
                body=[0]*16
                for j,i in enumerate(positions):body[i]=6*(1-2*((w>>j)&1))
                rows.append(body+[0])
    return rows


def main():
    parser=argparse.ArgumentParser(__doc__)
    parser.add_argument('--families',default='six,dual,flat,newshort')
    parser.add_argument('--blockers',type=int,default=8)
    parser.add_argument('--candidates',type=int,default=12000)
    parser.add_argument('--seconds',type=float,default=90)
    parser.add_argument('--pool-seed',type=int)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    if (args.out/'report.json').exists():raise ValueError('Use a new experiment directory')
    start=monotonic();source,group,generators,masks=seed_and_group()
    reps=representatives(group,masks,args.families.split(','))
    print('Seed',len(source),'representatives',len(reps),'group',len(group),flush=True)
    context=RunContext(args.out,{'max_materialized_points':max(20000,len(source)+args.candidates),
        'max_pairs':2000000000,'max_codewords':2000000,'max_rank_work':10000000000,
        'max_gram_entries':30000000,'solver_seconds':max(120,args.seconds)})
    artifact=Artifact('spherical_code',{'coordinates':source,'weights':[1]*16+[2],'norm_squared':288})
    seed_check=verify_artifact(artifact,context)
    assert seed_check.status=='exact'
    fingerprints=source_fingerprints()
    def progress(event):
        (args.out/'checkpoint.json').write_text(json.dumps(event,indent=2)+'\n')
        print(json.dumps({k:v for k,v in event.items() if k!='blocker_histogram_orbits'}),flush=True)
    result=select_partial_orbits(source,[1]*16+[2],288,reps,generators,context,
        blocker_limit=args.blockers,candidate_limit=args.candidates,seconds=args.seconds,
        pool_seed=args.pool_seed,on_progress=progress)
    report={'dimension':17,'baseline':5730,'result':result.facts['cardinality'],
        'families':args.families.split(','),'source_sha256':digest(artifact.data),
        'search':result.data['search'],'verification':result.verification.to_dict(),
        'engine_sources':fingerprints,
        'elapsed_seconds':monotonic()-start,'new_record':False,
        'scope':'Finite catalogue and explicit blocker/candidate limits; no exact upper bound'}
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    if result.facts['cardinality']>5730:
        (args.out/'candidate.json').write_text(json.dumps(result.to_dict())+'\n')
    print(json.dumps({k:report[k] for k in ('dimension','baseline','result','elapsed_seconds')}),flush=True)
    print(json.dumps(report['search']),flush=True)


if __name__=='__main__':main()
