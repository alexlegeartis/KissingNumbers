"""Exact theta, design-moment and finite-tail proof of the 77--79D bounds.

Literature inputs: Gamma72 is even unimodular of minimum 8, and every
nonempty shell of an extremal even unimodular 72D lattice is an 11-design.
No optimization output or floating-point angle comparison is trusted.
The large independent set is proved to exist, not supplied as a coordinate list.
"""
from collections import Counter
from fractions import Fraction as Q
from math import comb, prod
from pathlib import Path
import hashlib
import json
import numpy as np
import sympy as sp

HERE=Path(__file__).resolve().parent


def multiply(a,b):
    return [sum(a[j]*b[i-j] for j in range(i+1)) for i in range(6)]


def power(a,n):
    result=[1,0,0,0,0,0]
    for _ in range(n): result=multiply(result,a)
    return result


def theta_coefficients():
    # Theta has weight 36 and lies in span E4^9, E4^6 Delta,
    # E4^3 Delta^2, Delta^3. Minimum 8 sets coefficients q^1..q^3 to zero.
    e4=[1]+[240*sum(d**3 for d in range(1,n+1) if n%d==0) for n in range(1,6)]
    delta=[0,1,0,0,0,0]
    for n in range(1,6):
        factor=[0]*6
        for j in range(6//n+1):
            if n*j<6: factor[n*j]=(-1)**j*comb(24,j)
        delta=multiply(delta,factor)
    theta=power(e4,9)
    coefficients=[]
    for j in range(1,4):
        term=multiply(power(e4,9-3*j),power(delta,j))
        assert term[j]==1 and not any(term[:j])
        coefficient=-theta[j]
        coefficients.append(coefficient)
        theta=[x+coefficient*y for x,y in zip(theta,term)]
    assert theta[:4]==[1,0,0,0]
    return theta,coefficients


def solve_moments(population,radius_product,known):
    matrix=sp.Matrix([[j**(2*r) for j in range(6)] for r in range(6)])
    rhs=[]
    for r in range(6):
        moment=Q(population*prod(range(1,2*r,2))*radius_product**r,
                 prod(72+2*j for j in range(r)))
        rhs.append(sp.Rational(moment.numerator,moment.denominator)-
                   sum(c*j**(2*r) for j,c in known.items()))
    answer=matrix.inv()*sp.Matrix(rhs)
    assert all(v.is_Integer and v>=0 for v in answer)
    return {j:int(v) for j,v in enumerate(answer)}|known


def check_tails(z,rank,q):
    assert z.ndim==2 and z.shape[1]==8 and len(z)>0
    assert len(np.unique(z,axis=0))==len(z)
    rows={tuple(v) for v in z}
    assert all(tuple(-v) in rows for v in z)
    assert np.all(np.sum(z*z,axis=1)==q)
    assert sp.Matrix(z.tolist()).rank()==rank
    # Rank is checked exactly; all points lie in a Euclidean rank-k subspace.


def main():
    theta,coefficients=theta_coefficients()
    n8,n10=theta[4:6]
    assert n8==6218175600 and n10==15281788354560
    small=solve_moments(n8,80,{})
    # |y.u|<=5 for norm-ten y and minimal u. Norm-ten x with y.x=6
    # correspond bijectively to minimal u=y-x with y.u=4.
    shell=solve_moments(n10,100,{6:small[4],10:2})
    assert sum(shell.values())==n10
    assert all(c%2==0 for j,c in shell.items() if j)
    source=HERE/'norm10_tails.npz'
    rows={}
    with np.load(source) as data:
        for rank,q,scale in [(5,40,Q(16,25)),(6,40,Q(16,25)),(7,32,Q(12,19))]:
            z=data['tails%d'%rank].astype(np.int64)
            poles=data['poles%d'%rank].astype(np.int64)
            check_tails(z,rank,q)
            check_tails(poles,rank,8)
            assert sp.Matrix(np.vstack((z,poles)).tolist()).rank()==rank
            pd=poles@poles.T
            np.fill_diagonal(pd,0)
            assert int(pd.max())<=4
            assert scale<=Q(16,25) and scale>=Q(3,5)
            # Scale<=4/5 gives equatorial products<=4. Scale^2>=3/5
            # makes every pure-tail point compatible by Cauchy--Schwarz.
            numerator,denominator=scale.numerator,scale.denominator
            tail_coefficient=8*denominator-10*numerator
            profiles=Counter()
            for start in range(0,len(z),128):
                dots=z[start:start+128]@z.T
                counts=[np.count_nonzero(tail_coefficient*dots>(4*denominator-d*numerator)*q,axis=1)
                        for d in [10,4,5,6]]
                # For distinct head lines |D|<=3 there can be no conflict.
                assert 3*scale+(8-10*scale)<=4
                profiles.update(zip(*[v.tolist() for v in counts]))
            reciprocal=Q(0)
            for profile,count in profiles.items():
                self_closed,b4,b5,b6=profile
                closed=self_closed+sum(shell[d]//2*b for d,b in zip([4,5,6],[b4,b5,b6]))
                reciprocal+=Q(count,closed)
            alpha=(n10//2)*reciprocal
            pairs=-(-alpha.numerator//alpha.denominator)
            rows[str(72+rank)]={'bound':n8+2*pairs+len(poles),'equator':n8,
                'head_scale_squared':str(scale),'tail_squared_norm':q,'tail_vectors':len(z),
                'guaranteed_antipodal_pairs':pairs,'pure_tail_points':len(poles),
                'caro_wei_numerator':str(alpha.numerator),'caro_wei_denominator':str(alpha.denominator),
                'profiles':[list(p)+[c] for p,c in sorted(profiles.items())]}
    result={'status':'PASS','certificate_type':'Exact finite-graph existence proof; no full head list',
        'theta_coefficients':theta,'modular_form_coefficients':coefficients,
        'minimal_shell_absolute_products_against_norm10':small,
        'norm10_absolute_product_distribution':shell,'dimensions':rows,
        'tails_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'literature_inputs':['Gamma72 even unimodular, minimum 8 (Nebe 2010)',
                             'All shells are spherical 11-designs (Venkov; Nebe survey Theorem 4.1)',
                             'Classical modular-form description of even unimodular theta series']}
    (HERE/'norm10_bound_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__': main()
