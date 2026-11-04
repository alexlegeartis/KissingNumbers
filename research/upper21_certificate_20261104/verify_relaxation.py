"""Exact feasible pseudo-moments at 971, proving these relaxations stay loose.

This is not a point configuration. It certifies feasibility of the tested
two-point SDP with conditional maximal-clique and clique-product inequalities.
"""
import argparse,json,time
from fractions import Fraction as F
import numpy as np
from verify_upper import setup,D


def dot(c,x):return sum((int(a)*b for a,b in zip(c,x)),F(0))

def positive_definite(M):
    n=len(M);L=[[F(int(i==j)) for j in range(n)] for i in range(n)];diag=[]
    for i in range(n):
        di=M[i][i]-sum((L[i][j]**2*diag[j] for j in range(i)),F(0))
        assert di>0
        diag.append(di)
        for r in range(i+1,n):
            L[r][i]=(M[r][i]-sum((L[r][j]*L[i][j]*diag[j] for j in range(i)),F(0)))/di
    return min(diag)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--generate',action='store_true');a=p.parse_args()
    started=time.monotonic();O,P,z=setup();nv=len(z['reps']);target=971
    path=D/'relaxation971_rational.json'
    if a.generate:
        raw=np.load(D/'labels_cliques_four_products_target971_moments.npy')
        x=[F(round(float(t)*10**9),10**9) for t in raw]
        for j in np.flatnonzero(z['edges']):x[j]=F(0)
        x[int(O[0,0,0])]=F(1)
        row=np.bincount(O[:,0,:].ravel(),minlength=nv)
        free=[j for j in range(nv) if not z['edges'][j] and j!=O[0,0,0]]
        adjust=max(free,key=lambda j:row[j]);x[adjust]+=(target-dot(row,x))/int(row[adjust])
    else:x=list(map(F,json.loads(path.read_text())['moments']))
    assert len(x)==nv and min(x)>=0 and max(x)<=1
    assert x[O[0,0,0]]==1 and all(x[j]==0 for j in np.flatnonzero(z['edges']))
    assert dot(np.bincount(O[:,0,:].ravel(),minlength=nv),x)==target
    # Reconstruct each clique directly, rather than trusting saved coefficient rows.
    cliques=0;four=0
    for d in range(1024):
        base=sum((x[j] for j in O[d,0]),F(0))
        for p in range(21):
            value=base+x[O[d^P[p],0,p]]
            assert F(23)-F(21504,target)<=value<=1
            cliques+=1
            for q in range(21):
                if p==q:continue
                value=x[O[d,0,p]]+x[O[d,0,q]]+x[O[d^P[p],0,p]]+x[O[d^P[p],0,q]]
                assert value<=1;four+=1
    # All maximal-clique product types are covered by these verified pair orbits.
    for rep in z['reps']:
        d,p,q=map(int,(rep//441,(rep%441)//21,rep%21))
        base=sum((x[j] for j in O[d].ravel()),F(0))
        left=sum((x[j] for j in O[d^P[p],p]),F(0))
        right=sum((x[j] for j in O[d^P[q],:,q]),F(0))
        cross=x[O[d^P[p]^P[q],p,q]]
        for value,size in ((base,42),(base+left,43),(base+right,43),(base+left+right+cross,44)):
            assert F(size)-F(21504,target)<=value<=size//2
    # Check all 1024 characters, grouping equal coefficient blocks to avoid
    # repeating identical rational LDL factorizations.
    groups={}
    for chi in range(1024):
        signs=np.array([1-2*((d&chi).bit_count()%2) for d in range(1024)],dtype=np.int64)
        # We only need eight representatives because graph automorphisms permute
        # the characters. Check coverage under the independently verified group.
        groups[chi]=signs
    gens=np.load(D/'symmetry.npz')['generators'];covered=set();minpivot=None;blocks=0
    for chi in z['characters']:
        chi=int(chi);orbit={chi};frontier=[chi]
        for c in frontier:
            for g in gens:
                v=sum(((c&int(g[1<<i])).bit_count()%2)<<i for i in range(10))
                if v not in orbit:orbit.add(v);frontier.append(v)
        assert not covered.intersection(orbit);covered.update(orbit)
        B=np.array([np.bincount(O[:,a,b],weights=groups[chi],minlength=nv)
                    for a in range(21) for b in range(21)],dtype=np.int64).reshape(21,21,nv)
        M=[[dot(B[a,b],x) for b in range(21)] for a in range(21)]
        pivot=positive_definite(M);minpivot=pivot if minpivot is None else min(minpivot,pivot);blocks+=1
    assert len(covered)==1024
    if a.generate:path.write_text(json.dumps(dict(target=target,moments=list(map(str,x))),indent=2)+'\n')
    report=dict(status='PASS',relaxation_feasible_cardinality=target,actual_configuration=False,
                exact_psd_blocks=blocks,characters_covered=len(covered),cliques_checked=cliques,
                four_cliques_checked=four,minimum_ldl_pivot_approx=float(minpivot),
                arithmetic='Python Fraction; positive rational LDL pivots',seconds=time.monotonic()-started)
    (D/'relaxation_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
