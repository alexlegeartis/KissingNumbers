"""Exact replay of complete big families plus explicitly selected tiny points.

Run inside the portable results directory. NumPy is the only dependency.
The inherited population identities and published parent minimum/design
properties have the same scope as in verify_families.py.
"""
from pathlib import Path
import json,time,hashlib
import numpy as np
import verify_families as vf
ROOT=Path(__file__).resolve().parent

def main():
    start=time.monotonic();g,residues,mu,tiny=vf.geometry()
    m=dict(np.load(ROOT/'configuration_43.npz'));base=vf.configuration(ROOT/'configuration_43.npz',g,residues,mu,tiny)
    assert np.all(m['r3']>=12)
    raw=np.load(ROOT/'tiny_bodies.npz')['H'].astype(np.int64)
    bq=np.load(ROOT/'body_complement.npz');qi3=vf.inverse(bq['q'],3)
    co=vf.mm(vf.mm(raw,bq['b'].T),qi3);assert not np.any(co%288);co//=288
    lookup={v.tobytes():i for i,v in enumerate(residues)}
    labels=np.array([lookup[v.tobytes()] for v in co%3])
    by={int(c):np.flatnonzero(labels==c) for c in np.unique(labels)}
    a=np.load(ROOT/'explicit_points_43.npz');ix=a['body_indices'];z=a['head_numerators'].astype(object);D=int(a['head_denominator'])
    assert D>0 and ix.ndim==1 and len(ix)==len(z) and np.all((ix>=0)&(ix<len(raw)))
    P=m['span_projector_numerator'].astype(object);pd=int(m['span_projector_denominator']);assert np.array_equal(z@P,pd*z)
    X=raw[ix];c=labels[ix];C=m['classes'];R=m['r3'];Z=m['head_numerators'].astype(object);d=int(m['head_denominator'])
    N=10*D*D+np.sum(z*z,axis=1);M=R.astype(object)*d*d+np.sum(Z*Z,axis=1)
    B=np.empty((len(ix),len(C)),np.int64);refined=0
    for cl in np.unique(c):
        rows=np.flatnonzero(c==cl);dc=((g['digits'][C]-g['digits'][cl])%3)@g['powers'];md=mu[dc]
        B[rows]=10+R-18+md
        for j in np.flatnonzero(md==8):
            # All norm-10/3 differences in this class occur in the complete
            # explicit tiny pool. If no such difference yields the old norm,
            # the next possible squared difference norm is at least 16/3.
            contact=np.any(vf.mm(X[rows],raw[by[int(dc[j])]].T)==(int(R[j])-20)*144,axis=1)
            B[rows[~contact],j]-=6;refined+=int(np.sum(~contact))
    doubled=B.astype(object)*D*d+2*z@Z.T
    assert np.all((doubled<=0)|(doubled*doubled<=N[:,None]*M)), 'explicit/big angle'
    checks=len(ix)*len(C)
    for lo in range(0,len(ix),64):
        hi=min(lo+64,len(ix));body=vf.mm(X[lo:hi],X.T).astype(object)
        doubled=body*D*D+288*z[lo:hi]@z.T
        bad=(doubled>0)&(doubled*doubled>144**2*N[lo:hi,None]*N)
        bad[np.arange(hi-lo),np.arange(lo,hi)]=False
        assert not np.any(bad),('explicit/explicit angle',lo)
    checks+=len(ix)*(len(ix)-1)//2
    pure=m['pure_directions'].astype(object);pn=np.sum(pure*pure,axis=1);doubled=2*z@pure.T
    assert np.all((doubled<=0)|(doubled*doubled<=N[:,None]*pn)), 'explicit/pure angle'
    checks+=len(ix)*len(pure)
    result=dict(status='PASS',dimension=43,points=base['points']+len(ix),complete_families=len(C),explicit_tiny_points=len(ix),pure_points=len(pure),
                exact_angle_checks=base['exact_angle_checks']+checks,refined_body_bounds=refined,
                population_scope='Inherited population identities for complete families; explicit tiny vectors and every angle checked.',seconds=time.monotonic()-start)
    (ROOT/'replay.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))

if __name__=='__main__':main()
