#!/usr/bin/env python3
"""Read-only exact direction, pole and symbolic lift replay for dimensions 49–60.

SPDX-License-Identifier: MIT
Copyright (c) 2026 Alexey Kravatskiy

Input formats and root descriptions follow KissingNumbers' dim54-59-p48-coverage
verify.py and dim60-antipodal840-lift verify_tail.py/verify_poles.py. See
source/TAIL_SOURCE_PROVENANCE.json and source/LICENSE for attribution and hashes.
No network, optimizer, expanded P48 shell, or external filesystem input is used.
"""
from fractions import Fraction
from pathlib import Path
import argparse, hashlib, itertools, json, math
import numpy as np

HERE=Path(__file__).resolve().parent

def integer_array(a):
    assert a.dtype.kind in 'iu' and a.size
    assert int(a.min())>=-2**30 and int(a.max())<=2**30
    return a.astype(np.int64)

def abmax(a):
    return max(abs(int(a.min())),abs(int(a.max())))

def quadratic_sign(a,b,r):
    """Sign of a+b sqrt(r), entirely in arbitrary-precision integers."""
    a,b,r=int(a),int(b),int(r)
    if b==0:return (a>0)-(a<0)
    if a==0:return (b>0)-(b<0)
    if (a>0)==(b>0):return 1 if a>0 else -1
    difference=a*a-r*b*b
    return ((difference>0)-(difference<0))*(1 if a>0 else -1)

def check_code(A,B,c,r,N,n,k,relation=None):
    A,B,c=map(integer_array,(A,B,c));r,N=int(r),int(N)
    assert A.shape==B.shape and A.shape[0]==n and c.shape==(A.shape[1],)
    d=A.shape[1]
    assert N>0 and np.all(c>0) and r>=0
    assert (r==0 and not B.any()) or (r>0 and math.isqrt(r)**2!=r)
    if relation=='A2':
        assert d==3 and k==2 and np.all(A.sum(1)==0) and np.all(B.sum(1)==0)
    elif relation=='E7':
        assert d==8 and k==7 and np.all(A[:,6]+A[:,7]==0) and not B.any()
    elif relation=='E6':
        # Independent equations, with pivots in columns 6 and 5.
        assert d==8 and k==6 and np.all(A[:,6]+A[:,7]==0)
        assert np.all(A[:,5]+A[:,7]==0) and not B.any()
    else:assert d==k
    # Bound products, summation intermediates, scaled matrices and doubled
    # comparisons. Converting to Python int in sign avoids squaring overflow.
    bound=d*abmax(c)*(abmax(A)+abmax(B))**2*max(1,r)
    assert 4*bound<2**63 and 2*N<2**63
    U=(A*c)@A.T+r*((B*c)@B.T)
    V=(A*c)@B.T+(B*c)@A.T
    assert np.all(np.diag(U)==N) and not np.diag(V).any()
    keys={tuple(int(v) for v in np.r_[a,b]) for a,b in zip(A,B)}
    assert len(keys)==n and n%2==0
    assert all(tuple(-int(v) for v in np.r_[a,b]) in keys for a,b in zip(A,B))
    contacts=0
    for i in range(n):
        for j in range(i):
            s=quadratic_sign(2*int(U[i,j])-N,2*int(V[i,j]),r)
            assert s<=0,('tail violation',k,i,j)
            contacts+=s==0
    # Equal positive rays are impossible because fixed norm and pair <=1/2.
    return dict(tail_dimension=k,points=n,antipodal_pairs=n//2,pairs_checked=n*(n-1)//2,contacts=contacts,arithmetic='bounded int64 quadratic products and Python-integer sign decisions',embedding='weighted real coordinates sqrt(c)*(A+B*sqrt(r))/sqrt(N)',dimension_reason=relation or 'explicit coordinate count')

def simple_code(k):
    if k==1:return np.array([[-1],[1]],dtype=np.int64),1,None
    if k==2:
        rows=[]
        for i in range(3):
            for j in range(3):
                if i!=j:
                    v=np.zeros(3,dtype=np.int64);v[i]=1;v[j]=-1;rows.append(v)
        return np.array(rows),2,'A2'
    if k<=5:
        rows=[]
        for i,j in itertools.combinations(range(k),2):
            for s,t in itertools.product((-1,1),repeat=2):
                v=np.zeros(k,dtype=np.int64);v[i]=s;v[j]=t;rows.append(v)
        return np.array(rows),2,None
    rows=[]
    for i,j in itertools.combinations(range(8),2):
        for s,t in itertools.product((-2,2),repeat=2):
            v=np.zeros(8,dtype=np.int64);v[i]=s;v[j]=t;rows.append(v)
    rows.extend(np.array(v,dtype=np.int64) for v in itertools.product((-1,1),repeat=8) if sum(x<0 for x in v)%2==0)
    A=np.array(rows,dtype=np.int64)
    if k<8:A=A[A[:,6]+A[:,7]==0]
    if k<7:A=A[A[:,5]+A[:,7]==0]
    return A,8,{6:'E6',7:'E7',8:None}[k]

def rational_poles(path):
    raw=path.read_bytes()
    rows=[[Fraction(s) for s in line.split()] for line in raw.decode().splitlines() if line.strip()]
    assert len(rows)==841 and all(len(row)==12 for row in rows)
    # A separate common denominator for each ray preserves its unit direction.
    integers=[]
    for row in rows:
        scale=math.lcm(*(v.denominator for v in row))
        integers.append([v.numerator*(scale//v.denominator) for v in row])
    X=np.array(integers,dtype=object);norms=[sum(v*v for v in row) for row in integers]
    assert all(n>0 for n in norms)
    H=X@X.T
    checked=0
    for i in range(len(rows)):
        assert H[i,i]==norms[i]
        for j in range(i):
            dot=int(H[i,j])
            # For nonpositive dot the upper kissing bound is automatic.
            assert dot<=0 or 4*dot*dot<=norms[i]*norms[j],('pole violation',i,j)
            checked+=1
    return dict(points=841,dimension=12,pairs_checked=checked,normalization='Each exact printed rational ray divided by its own Euclidean norm',arithmetic='arbitrary-precision integers',distinctness='Every distinct input pair has normalized inner product <=1/2; hence no repeated positive ray',sha256=hashlib.sha256(raw).hexdigest())

def lift_proof():
    h=Fraction(1,2);a=Fraction(3,4);b=Fraction(1,4)
    values={'retained_equator_pair':h,'different_classes':a*h+b*h,'same_class_distinct_signed_heads':a*Fraction(1,3)+b,'same_signed_head_opposite_tails':a-b,'opposite_signed_heads':-a+b,'mixed_pole':h,'pole_pair':h}
    assert a+b==1 and all(v<=h for v in values.values())
    # Equator–mixed upper bound sqrt(3)/4 is positive and square3/16<1/4.
    assert a*h*h<h*h
    return dict(status='PASS',rational_pair_bounds={key:str(value) for key,value in values.items()},equator_mixed='sqrt(3)/4 < 1/2, checked by exact squaring',mixed_pole='Cauchy–Schwarz bounds tail half times any unit pole by1/2; pole code independent of direction code',ownership='Each selected signed mother head removed from equator and assigned to exactly one class; antipodal certificate uses both signs of each unique line',count='M - 2L + 4L + P = 52416000 + 2L + P',scope='Algebraic pair-category proof, not expansion or pairwise comparison of the 50-million-point configuration')

def verify(source=None):
    source=Path(source) if source is not None else HERE/'source'
    counts={1:2,2:6,3:12,4:24,5:40,6:72,7:126,8:240,9:306,10:510,11:604,12:840}
    reports=[];hashes={}
    for k,n in counts.items():
        if k<=8:
            A,N,relation=simple_code(k);B=np.zeros_like(A);c=np.ones(A.shape[1],dtype=np.int64);r=2
        else:
            name=f'capdirs_{k}.npz' if k<12 else 'tail840_exact.npz'
            path=source/name;raw=path.read_bytes();hashes[name]=hashlib.sha256(raw).hexdigest()
            # Read precisely the bytes whose hash is returned.
            import io
            z=np.load(io.BytesIO(raw),allow_pickle=False)
            A=z['A'];B=z['B'];relation=None
            if k<12:c=z['c'];r=int(z['R'].item());N=int(z['norm'].item())
            else:c=np.ones(12,dtype=np.int64);r=2;N=16
        report=check_code(A,B,c,r,N,n,k,relation)
        report['ambient_dimension']=48+k
        report['pole_points']=n if k<12 else 841
        report['pole_code']='same normalized tail code' if k<12 else 'independent rational-ray code'
        reports.append(report)
    poles=rational_poles(source/'coordinates841.txt')
    hashes['coordinates841.txt']=poles['sha256']
    return dict(status='PASS',tails=reports,independent_poles=poles,lift=lift_proof(),source_sha256=hashes,inherited_inputs='Identification of published integral P48 Gram lattice, minimum squared norm6 and shell population52416000. No mother-shell enumeration is performed.',head_scope='Final head ownership, norms and within-class pair inequalities require the separate head certificate replay.')

if __name__=='__main__':
    if not __debug__:raise SystemExit('Run without -O: assertions must be enabled.')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=HERE/'source')
    args=parser.parse_args()
    print(json.dumps(verify(args.source),indent=2))
