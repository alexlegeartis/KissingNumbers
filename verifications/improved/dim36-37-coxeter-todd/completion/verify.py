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
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "phase"))
from exact_lattice import enumerate_short, determinant, multiply, transpose

HERE=Path(__file__).resolve().parent


def main():
    started=time.monotonic()
    path=HERE/'certificate.json';cert=json.loads(path.read_text())
    assert cert['format']=='k12_phase_and_lift_orbits_v2'
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
    lifts=[tuple(row) for row in cert['lift_phase_orbits']]
    assert cert['lift_heads']==[2,-2] and len(set(lifts))==len(lifts)
    # Label 252 denotes the zero block. Fix the phase of the first nonzero
    # block to zero; the other nonzero block has one of three relative phases.
    lift_cells=set()
    for a,b,c,p,q in lifts:
        assert all(type(x) is int for x in (a,b,c,p,q))
        assert 0<=p<3 and 0<=q<3
        labels=(a,b,c); phases=(0,p,q)
        assert labels.count(252)==1 and all(0<=x<=252 for x in labels)
        nonzero=[i for i,x in enumerate(labels) if x!=252]
        assert phases[nonzero[0]]==0 and phases[labels.index(252)]==0
        lift_cells.add(tuple(0 if x==252 else roots[x] for x in labels))
    lift_base_pairs=0
    for a,b,c in lift_cells:
        for x,y,z in base:
            assert D[a][x]+D[b][y]+D[c][z]>=8
            lift_base_pairs+=1
    # A weight-four block has dot <=2 against a block with a different
    # nonzero coset label. Thus a lift/equator dot can exceed 6 only when
    # both of the lift's nonzero labels coincide with the equator labels.
    T.append([[0]*3 for _ in range(252)])
    for row in T:row.append([0]*3)
    lift_equator_pairs=0
    for u in lifts:
        nonzero=[i for i in range(3) if u[i]!=252]
        i,j=nonzero
        candidates=set(rows[i][u[i]]).intersection(rows[j][u[j]])
        for index in candidates:
            v=orbits[index]
            maximum=max(T[u[0]][v[0]][t]+T[u[1]][v[1]][(v[3]-u[3]+t)%3]
                        +T[u[2]][v[2]][(v[4]-u[4]+t)%3] for t in range(3))
            assert maximum<=6
            lift_equator_pairs+=1
    lift_pairs=0
    for index,u in enumerate(lifts):
        for v in lifts[:index]:
            maximum=max(T[u[0]][v[0]][t]+T[u[1]][v[1]][(v[3]-u[3]+t)%3]
                        +T[u[2]][v[2]][(v[4]-u[4]+t)%3] for t in range(3))
            assert maximum<=2
            lift_pairs+=1
    # Body norm 8 plus head norm 4 gives norm 12. Same-height body
    # products <=2 give total products <=6; opposite heads give <=8-4=4.
    # Each orbit has three different bodies, so the two heights add six points.
    n37=n36+6*len(lifts)
    assert n37==cert['points_37']
    result=dict(points_36=n36,points_37=n37,base_points=base_points,
                phase_orbits=len(orbits),phase_cells_used=len(cells),
                lift_phase_orbits=len(lifts),lift_cells_used=len(lift_cells),
                base_coset_pairs_checked=base_pairs,
                phase_base_coset_pairs_checked=core_pairs,
                phase_orbit_pairs_checked=phase_pairs,
                lift_base_coset_pairs_checked=lift_base_pairs,
                lift_equator_orbit_pairs_checked=lift_equator_pairs,
                lift_orbit_pairs_checked=lift_pairs,
                exact_shell_counts=dict(shell_counts),squared_norm=12,
                minimum_squared_distance_at_least=12,
                unit_sphere_max_inner_product='1/2',all_checks_pass=True,
                certificate_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                arithmetic='Python integers and fractions.Fraction',
                elapsed_seconds=time.monotonic()-started)
    (HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
