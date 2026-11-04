"""Rational dual certificates for freezing half of a recovered kernel.

Average an arbitrary dominating set over a 32-word subgroup T of the kernel.
If T is entirely omitted, the quotient indicator at zero is one. Nonnegative
weights y on as-yet-undominated quotient vertices, with sum(y)=1 and all
other-column loads <=1, force another unit of quotient hole mass: >=64 holes.
This does not assume periodicity of the new configuration.
"""
import argparse,json
from fractions import Fraction
import numpy as np
from model import C,N,DATA
from quotients import span_basis


def quotient(mask):
    H=np.load(DATA/'recovery.npz')['holes'][0]
    K,_=span_basis(C[H]); lookup={int(w):i for i,w in enumerate(C)}
    T=K[[not ((i&mask).bit_count()%2) for i in range(64)]]
    reps,orb=np.unique(np.min(C[:,None]^T,axis=1),return_inverse=True)
    QN=orb[N[[lookup[int(w)] for w in reps]]]
    A=np.eye(32,dtype=int)
    for u in range(32):
        for v in QN[u]: A[u,v]+=1
    assert len(T)==32 and len(reps)==32 and A[0,0]==1
    return A,T


def verify(records):
    for r in records:
        A,T=quotient(r['mask'])
        y=[Fraction(s) for s in r['weights']]
        assert len(y)==32 and min(y)>=0 and sum(y)==1
        assert all(not y[i] or A[i,0]==0 for i in range(32))
        assert all(sum(y[i]*int(A[i,j]) for i in range(32))<=1 for j in range(1,32))
        assert sorted(map(int,T))==r['fixed_subgroup']
    return dict(status='PASS',rational_certificates=len(records),
                translated_fixed_half_sets=2*len(records),
                translated_cosets_in_full_component=32*len(records),
                conditional_hole_lower_bound=64,
                scope='Any dominating set containing any of the certified fixed 32-hole halves')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify-only',action='store_true');a=p.parse_args()
    path=DATA/'half_certificates.json'
    if a.verify_only:
        records=json.loads(path.read_text())['certificates']
    else:
        from scipy.optimize import linprog
        records=[];failed=[]
        for mask in range(1,64):
            A,T=quotient(mask); rows=np.flatnonzero(A[:,0]==0)
            lp=linprog(np.zeros(len(rows)),A_ub=A[rows,1:].T,b_ub=np.ones(31),
                       A_eq=np.ones((1,len(rows))),b_eq=[1.],bounds=(0,None),method='highs')
            if not lp.success: failed.append(mask);continue
            y=[Fraction(0) for _ in range(32)]
            for i,x in zip(rows,lp.x): y[i]=Fraction(float(x)).limit_denominator(1000000)
            records.append(dict(mask=mask,fixed_subgroup=sorted(map(int,T)),weights=list(map(str,y))))
        verify(records)
        path.write_text(json.dumps(dict(certificates=records,uncertified_masks=failed),indent=2)+'\n')
    report=verify(records)
    (DATA/'halves_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)
