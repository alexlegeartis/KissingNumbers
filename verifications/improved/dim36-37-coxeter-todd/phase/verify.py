"""Independent exact verification of the 36D phase code and 37D extension.

Python 3.9+, standard library only. No search or solver is run.
The coordinate convention is row coefficients, norm x G x^T.
"""
from collections import Counter
from pathlib import Path
import hashlib
import json
import math
import time
from exact_lattice import enumerate_short, determinant, multiply, transpose

HERE=Path(__file__).resolve().parent


def main():
    started=time.monotonic()
    path=HERE/'certificate.json';cert=json.loads(path.read_text())
    assert cert['format']=='k12_phase_orbits_and_latitudes_v1'
    assert cert['squared_norm']==12
    G,M=cert['gram_12'],cert['sublattice_matrix']
    assert len(G)==len(M)==12 and all(len(r)==12 for r in G+M)
    assert all(isinstance(x,int) for r in G+M for x in r)
    assert G==transpose(G)
    assert multiply(multiply(M,G),transpose(M))==[[3*x for x in r] for r in G]
    assert abs(determinant(M))==729
    first,second=cert['label_first'],cert['label_second']
    assert len(first)==len(second)==6 and sorted(first+second)==list(range(12))
    powers=[3**i for i in range(6)]
    def label(v):return sum(((v[a]+v[b])%3)*p for a,b,p in zip(first,second,powers))
    assert all(label(r)==0 for r in M)
    # The kernel of label() and the row lattice of M both have index 729.
    short=list(enumerate_short(G))
    shell_counts=Counter(n for _,n in short)
    assert shell_counts=={4:756,6:4032,8:20412}
    weights=[100]*729;mult=[0]*729;weights[0],mult[0]=0,1
    minimum_vectors=[set() for _ in range(729)];minimum_vectors[0].add((0,)*12)
    for v,n in short:
        assert sum(v[i]*G[i][j]*v[j] for i in range(12) for j in range(12))==n
        c=label(v)
        if n<weights[c]:weights[c],mult[c],minimum_vectors[c]=n,1,{v}
        elif n==weights[c]:mult[c]+=1;minimum_vectors[c].add(v)
    assert Counter(zip(weights,mult))=={(0,1):1,(4,3):252,(6,18):224,(8,81):252}
    digits=[[i//p%3 for p in powers] for i in range(729)]
    D=[[weights[sum(((a-b)%3)*p for a,b,p in zip(x,y,powers))] for y in digits] for x in digits]
    base=[tuple(c) for c in cert['base_outer_cosets']]
    assert len(base)==1429 and base[0]==(0,0,0) and len(set(base))==len(base)
    assert all(len(c)==3 and all(isinstance(v,int) and 0<=v<729 for v in c) for c in base)
    assert all(sum(weights[v] for v in c)==12 and 0 in c for c in base[1:])
    base_pairs=0
    for i,(a,b,c) in enumerate(base):
        for x,y,z in base[:i]:
            assert D[a][x]+D[b][y]+D[c][z]>=12
            base_pairs+=1
    base_points=3*756+sum(math.prod(mult[x] for x in c) for c in base[1:])
    assert base_points==403704
    phi=[[int(i==j)-M[i][j] for j in range(12)] for i in range(12)]
    assert multiply(multiply(phi,G),transpose(phi))==G
    assert multiply(multiply(phi,phi),phi)==[[int(i==j) for j in range(12)] for i in range(12)]
    roots=[i for i,w in enumerate(weights) if w==4]
    zero=cert['root_phase_zero_vectors']
    assert len(zero)==252
    X=[]
    for i,v in enumerate(zero):
        assert len(v)==12 and all(isinstance(a,int) for a in v)
        x=multiply([v],phi)[0];y=multiply([x],phi)[0]
        phases=[tuple(v),tuple(x),tuple(y)]
        assert set(phases)==minimum_vectors[roots[i]]
        X.append(phases)
    weighted=multiply(zero,G)
    T=[[[sum(a*b for a,b in zip(weighted[i],X[j][t])) for t in range(3)]
        for j in range(252)] for i in range(252)]
    orbits=[tuple(c) for c in cert['phase_orbits']]
    assert len(set(orbits))==len(orbits)
    assert all(len(c)==5 and all(isinstance(x,int) for x in c)
               and all(0<=x<252 for x in c[:3]) and all(0<=x<3 for x in c[3:]) for c in orbits)
    cells=set(c[:3] for c in orbits)
    core_pairs=0
    for a,b,c in cells:
        da,db,dc=D[roots[a]],D[roots[b]],D[roots[c]]
        for x,y,z in base:
            assert da[x]+db[y]+dc[z]>=12
            core_pairs+=1
    # Different cells with three unequal labels have distance >= 4+4+4.
    # For every remaining pair, evaluate all three relative phase shifts.
    rows=[[[] for _ in range(252)] for _ in range(3)]
    phase_pairs=0;maximum_inner=-12
    for i,v in enumerate(orbits):
        a,b,c,p,q=v
        for coordinate in range(3):
            for index in rows[coordinate][v[coordinate]]:
                u=orbits[index]
                if coordinate>0 and a==u[0]:continue
                if coordinate>1 and b==u[1]:continue
                aa,bb,cc,pp,qq=u
                dp,dq=(pp-p)%3,(qq-q)%3
                value=max(T[a][aa][t]+T[b][bb][(dp+t)%3]+T[c][cc][(dq+t)%3] for t in range(3))
                assert value<=6,('phase conflict',i,index,value)
                maximum_inner=max(maximum_inner,value);phase_pairs+=1
            rows[coordinate][v[coordinate]].append(i)
    n36=base_points+3*len(orbits)
    assert n36==cert['points_36']
    # Each orbit has three distinct points. In its first block their phases
    # are 0,1,2, and its other phases are t+p,t+q. The tuple representation
    # is unique, so there is no duplicate between different orbits.
    lifts=[tuple(c) for c in cert['lift_body_cosets']]
    assert len(set(lifts))==len(lifts)
    assert cert['lift_heads']==[2,-2]
    assert all(len(c)==3 and all(isinstance(x,int) and 0<=x<729 for x in c)
               and sorted(weights[x] for x in c)==[0,4,4] for c in lifts)
    all_equator_cosets=base+[tuple(roots[x] for x in c) for c in cells]
    lift_pairs=0;lift_equator_pairs=0
    for i,(a,b,c) in enumerate(lifts):
        da,db,dc=D[a],D[b],D[c]
        for x,y,z in lifts[:i]:
            assert da[x]+db[y]+dc[z]>=12
            lift_pairs+=1
        for x,y,z in all_equator_cosets:
            assert da[x]+db[y]+dc[z]>=8
            lift_equator_pairs+=1
    # A lifted body has squared norm 8 and head +/-2, hence norm 12.
    # Same-head distances are >=12; opposite heads alone contribute 16.
    # For any equatorial norm-12 vector, body distance >=8 implies inner
    # product <= (8+12-8)/2 = 6; its head is zero.
    n37=n36+2*sum(math.prod(mult[x] for x in c) for c in lifts)
    assert n37==cert['points_37']
    result={'points_36':n36,'points_37':n37,'squared_norm':12,
            'minimum_squared_distance_at_least':12,'unit_sphere_max_inner_product':'1/2',
            'base_points':base_points,'phase_orbits':len(orbits),'phase_cells_used':len(cells),
            'lift_body_cells':len(lifts),'base_coset_pairs_checked':base_pairs,
            'phase_base_coset_pairs_checked':core_pairs,'phase_orbit_pairs_checked':phase_pairs,
            'maximum_inner_product_in_phase_check':maximum_inner,
            'lift_coset_pairs_checked':lift_pairs,'lift_equator_coset_pairs_checked':lift_equator_pairs,
            'exact_shell_counts':dict(shell_counts),'all_checks_pass':True,
            'certificate_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'arithmetic':'Python integers and fractions.Fraction; standard library only',
            'elapsed_seconds':time.monotonic()-started}
    (HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
