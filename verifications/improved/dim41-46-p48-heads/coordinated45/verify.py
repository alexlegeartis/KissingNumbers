"""Exact integer replay. Requires Python 3 and NumPy; no optimizer."""
from pathlib import Path
import json,time,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent
def verify():
    start=time.monotonic();m=np.load(ROOT/'configuration.npz');g=np.load(ROOT/'geometry.npz')
    import verify_inherited
    source,residues,source_mu,source_tiny=verify_inherited.geometry()
    assert np.array_equal(g['digits'],source['digits']) and np.array_equal(g['powers'],source['powers'])
    assert np.array_equal(g['mu'],source_mu) and np.array_equal(g['tiny'],source_tiny)
    assert np.array_equal(g['body_residues']%3,residues)
    assert np.array_equal(g['body_gram'],np.load(ROOT/'body_complement.npz')['q'])
    c=m['classes'];r=m['r3'];w=m['weights'];z=m['head_numerators'].astype(object);D=int(m['head_denominator']);pure=m['pure_directions'].astype(object)
    DIG=g['digits'];PW=g['powers'];MU=g['mu'];tiny=g['tiny']
    assert D>0 and len(c)==len(r)==len(w)==len(z)
    assert np.all(np.isin(r,[10,12,14,16,18])) and np.all(MU[c]%6==(18-r)%6)
    assert np.all(c[r==18]!=0) and np.all(MU[c[r==14]]==4) and np.all(MU[c[r==12]]==6)
    expected=np.where(r==18,1084590,np.where(r==16,116235-40*tiny[c],np.where(r==14,9180,np.where(r==12,495,tiny[c]))))
    assert np.array_equal(w,expected)
    A=m['span_projector_numerator'].astype(object);a=int(m['span_projector_denominator'])
    assert a>0 and np.array_equal(A,A.T) and np.array_equal(A@A,a*A) and np.trace(A)==5*a
    assert np.array_equal(z@A,a*z) and np.array_equal(pure@A,a*pure)
    N=r.astype(object)*D*D+np.sum(z*z,axis=1);assert np.all(N>0)
    checks=equalities=0
    for i in range(len(c)):
        dc=((DIG[c[i]]-DIG[c[:i+1]])%3)@PW
        b=(r[i]+r[:i+1]-18+MU[dc]).astype(object)
        b[(c[:i+1]==c[i])&(r[:i+1]==r[i])]=2*int(r[i]);b[-1]=2*int(r[i])-18
        dot=b*D*D+2*z[:i+1]@z[i];slack=N[i]*N[:i+1]-dot*dot
        assert not np.any((dot>0)&(slack<0)),('family',i)
        equalities+=int(np.sum((dot>0)&(slack==0)));checks+=i+1
    pn=np.sum(pure*pure,axis=1);assert np.all(pn>0)
    for i,v in enumerate(pure):
        dot=2*z@v;assert np.all((dot<=0)|(dot*dot<=pn[i]*N)),('pure/family',i)
        dot=2*pure[:i]@v;assert np.all((dot<=0)|(dot*dot<=pn[i]*pn[:i])),('pure/pure',i)
        checks+=len(c)+i
    added=np.load(ROOT/'added_bodies.npz');body_counts={}
    q=g['body_gram'];res=g['body_residues']
    for key in added.files:
        _,cl,rl=key.split('_');cl=int(cl);rl=int(rl);co=added[key].astype(np.int64)
        assert len(np.unique(co,axis=0))==len(co)
        assert np.all(co%3==res[cl]%3)
        assert np.all(np.sum((co@q)*co,axis=1)==3*rl)
        dots=co@q@co.T;np.fill_diagonal(dots,-100000)
        assert np.max(dots)<=3*rl-27
        assert len(co)==495 and rl==12
        body_counts[str(cl)]=len(co)
    result=dict(status='PASS',dimension=45,points=int(w.sum())+len(pure),families=len(c),pure_points=len(pure),exact_angle_checks=checks,equalities=equalities,explicit_added_body_counts=body_counts,seconds=time.monotonic()-start,scope='Exact head-angle and new-body replay; inherited K40 shell populations and coset distance bounds',sha256=hashlib.sha256((ROOT/'configuration.npz').read_bytes()).hexdigest())
    print(json.dumps(result,indent=2));return result
if __name__=='__main__':verify()
