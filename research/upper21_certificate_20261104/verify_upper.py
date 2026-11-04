"""Standalone verifier: NumPy and Python exact fractions; no optimization solver.

Checks exact geometry and automorphisms, independently checks connectivity
of every claimed moment orbit, and replays all rational proof terms.
"""
import json,time
from fractions import Fraction
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent;D=HERE/'data';OLD=HERE.parent/'compute19_21_20261104'/'data'

def setup():
    graph=np.load(OLD/'component21.npz');C=graph['component'];N=graph['neighbors'];P=N[0]
    assert np.array_equal(N,np.arange(1024)[:,None]^P)
    assert np.array_equal(C[np.arange(1024)[:,None]^np.arange(1024)],C[:,None]^C)
    geom=np.load(OLD/'baseline21.npz');profiles=geom['profiles'].astype(np.int64)
    signs=1-2*geom['code'][C,:21].astype(np.int64)
    assert np.all((profiles*profiles).sum(1)==164)
    dots=np.einsum('di,ai,bi->dab',signs,profiles,profiles,optimize=True)
    edge=dots>82
    for a in range(21):edge[0,a,a]=False
    z=np.load(D/'label_moments.npz');O=z['orbits'];nv=len(z['reps'])
    assert np.array_equal(z['edges'][O],edge)
    gens=np.load(D/'symmetry.npz')['generators']
    d=np.repeat(np.arange(1024),441);a=np.tile(np.repeat(np.arange(21),21),1024);b=np.tile(np.arange(21),1024*21)
    state_perms=[];label_perms=[]
    for g in gens:
        assert g[0]==0 and len(np.unique(g))==1024
        assert np.array_equal(np.sort(g[N],axis=1),np.sort(N[g],axis=1))
        for bit in (1<<np.arange(10)):assert np.array_equal(g[np.arange(1024)^bit],g^g[bit])
        lp=np.array([int(np.flatnonzero(P==g[p])[0]) for p in P]);label_perms.append(lp)
        perm=(g[d]*21+lp[a])*21+lp[b];state_perms.append(perm)
    state_perms.append((d*21+b)*21+a)
    perms=np.array(state_perms,dtype=np.int32);flat=O.ravel()
    for perm in perms:assert np.array_equal(flat,flat[perm])
    # Connectivity establishes that an orbit class has not merged inequivalent
    # pairs. Invariance alone would not establish this.
    seen=np.zeros(len(flat),bool)
    for j,rep in enumerate(z['reps']):
        assert flat[rep]==j
        frontier=np.array([rep],dtype=np.int32);seen[rep]=True;count=1
        while len(frontier):
            targets=np.unique(perms[:,frontier]);frontier=targets[~seen[targets]]
            seen[frontier]=True;count+=len(frontier)
        assert count==np.sum(flat==j)
    assert seen.all()
    labels={0}
    while True:
        new=labels|{int(g[p]) for g in label_perms for p in labels}
        if new==labels:break
        labels=new
    assert len(labels)==21
    return O,P,z


def verify(O,P,z,certificate):
    nv=len(z['reps']);acc=[Fraction(0)]*nv;rhs=Fraction(0);blocks={}
    def add(c,b,w):
        nonlocal rhs
        for j,value in enumerate(c):
            if value:acc[j]+=w*int(value)
        rhs+=w*b
    for term in certificate['inequalities']:
        w=Fraction(term['multiplier']);assert w>=0
        c=np.zeros(nv,dtype=np.int64);b=0
        if term['kind']=='nonnegative':c[term['variable']]=1
        elif term['kind']=='at_most_one':c[term['variable']]=-1;b=-1
        elif term['kind']=='clique':
            d,p=term['word'],term['line'];assert 0<=d<1024 and 0<=p<21
            c=-np.bincount(O[d,0],minlength=nv);c[O[d^P[p],0,p]]-=1;b=-1
        elif term['kind']=='psd':
            chi=term['character'];assert 0<=chi<1024
            if chi not in blocks:
                signs=np.array([1-2*((d&chi).bit_count()%2) for d in range(1024)],dtype=np.int64)
                blocks[chi]=np.array([np.bincount(O[:,a,b],weights=signs,minlength=nv)
                                     for a in range(21) for b in range(21)],dtype=np.int64).reshape(21,21,nv)
            v=np.array(term['vector'],dtype=np.int64);assert v.shape==(21,) and np.max(abs(v))<=4096
            c=np.einsum('a,abi,b->i',v,blocks[chi],v)
            div=term['divisor'];assert div>0 and np.all(c%div==0);c//=div
        else:raise AssertionError('unknown proof term')
        add(c,b,w)
    for term in certificate['equalities']:
        c=np.zeros(nv,dtype=np.int64);b=0
        if term['kind']=='diagonal':c[O[0,0,0]]=1;b=1
        elif term['kind']=='edge':
            j=term['variable'];assert z['edges'][j];c[j]=1
        elif term['kind']=='cardinality':
            c=np.bincount(O[:,0,:].ravel(),minlength=nv);b=term['value'];assert b==certificate['target']
        else:raise AssertionError('unknown equality')
        add(c,b,Fraction(term['multiplier']))
    assert all(a==0 for a in acc)
    assert rhs==Fraction(certificate['contradiction']) and rhs>0
    assert certificate['upper_component']==certificate['target']-1
    return rhs


if __name__=='__main__':
    start=time.monotonic();O,P,z=setup();c=json.loads((D/'upper971_certificate.json').read_text())
    gap=verify(O,P,z,c)
    # A changed cardinality equality must not replay as the same certificate.
    import copy
    bad=copy.deepcopy(c)
    for term in bad['equalities']:
        if term['kind']=='cardinality':term['value']-=1
    try:verify(O,P,z,bad)
    except AssertionError:control='PASS'
    else:raise AssertionError('corrupted certificate accepted')
    report=dict(status='PASS',excluded_size=c['target'],upper_component=c['upper_component'],
                family_total_upper=c['family_total_upper'],moment_orbits=len(z['reps']),
                inequalities=len(c['inequalities']),equalities=len(c['equalities']),
                exact_contradiction=str(gap),corruption_control=control,
                geometry_cases=1024*21*21,seconds=time.monotonic()-start,
                scope=c['scope'])
    (D/'upper_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)
