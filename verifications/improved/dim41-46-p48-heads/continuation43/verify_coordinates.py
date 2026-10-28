"""Independent coefficient-coordinate replay of the 43D tiny-layer continuation.

Uses integer products in the rank-40 body metric and integer square-root
cutoffs. The companion replay uses 48 ambient body coordinates and squaring.
Complete-family populations retain the stated P48n minimum/design inputs.
"""
import argparse
import hashlib
import json
from math import isqrt
from pathlib import Path

import numpy as np


def mm(a,b):
    a=np.asarray(a,dtype=np.int64);b=np.asarray(b,dtype=np.int64)
    assert a.shape[-1]*int(np.abs(a).max(initial=0))*int(np.abs(b).max(initial=0))<2**62
    return a@b


def inverse_mod3(a):
    n=len(a);v=np.c_[a.copy()%3,np.eye(n,dtype=np.int64)]
    for k in range(n):
        pivot=next(i for i in range(k,n) if v[i,k])
        v[[k,pivot]]=v[[pivot,k]];v[k]=(v[k]*int(v[k,k]))%3
        for i in range(n):
            if i!=k:v[i]=(v[i]-v[i,k]*v[k])%3
    assert np.array_equal(v[:,:n],np.eye(n,dtype=np.int64))
    return v[:,n:]


def accept(p,n,m):
    # Retain the sign before applying the exact integer square root.
    return int(p)<=0 or int(p)<=isqrt(int(n)*int(m))


def check(root,points,previous,delta):
    source=root/'partial43';model=dict(np.load(source/'configuration_43.npz',allow_pickle=False))
    geom=dict(np.load(source/'class_geometry.npz',allow_pickle=False))
    body=dict(np.load(source/'body_complement.npz',allow_pickle=False))
    raw=np.load(source/'tiny_bodies.npz',allow_pickle=False)['H'].astype(np.int64)
    q,b=body['q'].astype(np.int64),body['b'].astype(np.int64)
    assert q.shape==(40,40) and np.array_equal(q,q.T)
    assert np.array_equal(mm(b,b.T),96*q)
    # Floating inversion only proposes J; the integer identity certifies it.
    J=np.rint(3*np.linalg.inv(q.astype(float))).astype(np.int64)
    assert np.array_equal(mm(q,J),3*np.eye(40,dtype=np.int64))
    numerator=mm(mm(raw,b.T),J);assert not np.any(numerator%288)
    coefficients=numerator//288
    assert np.array_equal(mm(coefficients,b),raw)
    assert len(raw)==30240 and len(np.unique(coefficients,axis=0))==30240
    aq=mm(coefficients,q)
    assert 40*int(np.abs(aq).max())*int(np.abs(coefficients).max())<2**62
    assert np.all(np.sum(aq*coefficients,axis=1)==30)

    # Recover the linear gluing map from eight independently chosen family
    # residues, rather than importing the ambient-coordinate verifier.
    C=model['classes'];R=model['r3'];digits=geom['digits'];powers=geom['powers']
    pivots={};selected=[]
    for i,c in enumerate(C):
        row=digits[c].copy()%3
        for k,v in sorted(pivots.items()):row=(row-row[k]*v)%3
        nz=np.flatnonzero(row)
        if len(nz):
            k=int(nz[0]);pivots[k]=(row*int(row[k]))%3;selected.append(i)
        if len(selected)==8:break
    assert len(selected)==8
    inverse=inverse_mod3(digits[C[selected]])
    residues=((digits@inverse)%3)@model['family_residues'][selected]%3
    assert len(np.unique(residues,axis=0))==6561
    assert np.array_equal(residues[C],model['family_residues']%3)
    lookup={row.tobytes():i for i,row in enumerate(residues)}
    labels=np.array([lookup[row.tobytes()] for row in coefficients%3])
    tiny=np.bincount(labels,minlength=6561)
    assert np.array_equal(tiny,np.load(source/'population_table.npz')['tiny'])
    by={int(c):np.flatnonzero(labels==c) for c in np.unique(labels)}
    u=geom['u'].astype(np.int64);norm=np.sum(u*u,axis=1)
    assert np.all(u%2==u[:,:1]%2) and np.all(u.sum(axis=1)%4==0)
    assert np.array_equal(mm(geom['A'],geom['ai']),4*np.eye(8,dtype=np.int64))
    ui=mm(u,geom['ai']);assert not np.any(ui%4)
    ids=(((ui//4)%3)@powers)
    assert np.array_equal(ids,geom['ids']) and np.all(norm<=32)
    mu=np.full(6561,10000,dtype=np.int64);np.minimum.at(mu,ids,norm)
    assert np.array_equal(mu,geom['mu']) and not np.any(mu%4)
    mu//=4
    expected=np.where(R==18,np.where(C==0,1092000,1084590),np.where(
        R==16,116235-40*tiny[C],np.where(R==14,9180,495)))
    assert np.all(np.isin(R,[12,14,16,18])) and np.array_equal(expected,model['weights'])

    z=model['head_numerators'].astype(object);d=int(model['head_denominator'])
    family_n=R.astype(object)*d*d+np.sum(z*z,axis=1)
    assert d>0 and np.all(family_n>0)
    pure=model['pure_directions'].astype(object);pure_n=np.sum(pure*pure,axis=1)
    assert np.all(pure_n>0)
    assert model['dimension']==43 and z.shape==(len(C),3)
    assert np.array_equal(model['span_projector_numerator'],np.eye(3,dtype=np.int64))
    assert int(model['span_projector_denominator'])==1
    # Check the unchanged complete families and poles by a second arithmetic path.
    family_checks=0
    for i in range(len(C)):
        dc=((digits[C[:i+1]]-digits[C[i]])%3)@powers
        beta=R[:i+1]+R[i]-18+mu[dc]
        beta[(C[:i+1]==C[i])&(R[:i+1]==R[i])]=2*int(R[i]);beta[-1]=2*int(R[i])-18
        values=beta.astype(object)*d*d+2*z[:i+1]@z[i]
        assert all(accept(v,family_n[i],n) for v,n in zip(values,family_n[:i+1]))
        family_checks+=i+1
    for i,v in enumerate(pure):
        assert all(accept(a,pure_n[i],n) for a,n in zip(2*z@v,family_n))
        assert all(accept(a,pure_n[i],n) for a,n in zip(2*pure[:i]@v,pure_n[:i]))
        family_checks+=len(C)+i

    final=dict(np.load(points,allow_pickle=False));ix=final['body_indices'];h=final['head_numerators'].astype(object);D=int(final['head_denominator'])
    assert ix.dtype.kind in 'iu' and final['head_numerators'].dtype.kind in 'iu'
    assert ix.ndim==1 and h.shape==(len(ix),3) and D>0
    assert np.all((ix>=0)&(ix<len(raw)))
    assert len({(int(i),*map(int,row)) for i,row in zip(ix,h)})==len(ix)
    before=dict(np.load(previous,allow_pickle=False));change=dict(np.load(delta,allow_pickle=False))
    deleted=change['deleted_indices'];added=change['added_body_indices'];added_h=change['added_head_numerators']
    assert int(change['head_denominator'])==D and D%int(before['head_denominator'])==0
    assert np.array_equal(deleted,np.unique(deleted)) and np.all((deleted>=0)&(deleted<len(before['body_indices'])))
    kept=np.setdiff1d(np.arange(len(before['body_indices'])),deleted)
    restored_ix=np.r_[before['body_indices'][kept],added]
    restored_h=np.r_[before['head_numerators'][kept].astype(object)*(D//int(before['head_denominator'])),added_h]
    assert np.array_equal(ix,restored_ix) and np.array_equal(h,restored_h)
    N=10*D*D+np.sum(h*h,axis=1);a=coefficients[ix];ca=labels[ix]
    bounds=np.empty((len(ix),len(C)),dtype=np.int64);refined=0
    for c in np.unique(ca):
        rows=np.flatnonzero(ca==c);dc=((digits[C]-digits[c])%3)@powers
        bounds[rows]=10+R-18+mu[dc]
        for j in np.flatnonzero(mu[dc]==8):
            # A candidate difference t has (a+t)Q(a+t)=3*r_j exactly when
            # 2*aQt=3*(r_j-20). Exhaust all shortest differences in its class.
            products=2*mm(aq[ix[rows]],coefficients[by[int(dc[j])]].T)
            contact=np.any(products==3*(int(R[j])-20),axis=1)
            bounds[rows[~contact],j]-=6;refined+=int(np.sum(~contact))
    for lo in range(0,len(ix),64):
        hi=min(lo+64,len(ix));values=bounds[lo:hi].astype(object)*D*d+2*h[lo:hi]@z.T
        for k,row in enumerate(values):assert all(accept(v,N[lo+k],n) for v,n in zip(row,family_n))
    tiny_pairs=0
    for lo in range(0,len(ix),32):
        hi=min(lo+32,len(ix));dots=mm(aq[ix[lo:hi]],a[:hi].T)
        values=2*D*D*dots.astype(object)+6*h[lo:hi]@h[:hi].T
        for k,row in enumerate(values):
            i=lo+k
            assert all(int(v)<=0 or int(v)<=isqrt(9*int(N[i])*int(n)) for v,n in zip(row[:i],N[:i])),('tiny pair',i)
            tiny_pairs+=i
    for i,row in enumerate(h):assert all(accept(v,N[i],n) for v,n in zip(2*pure@row,pure_n))
    total=int(sum(map(int,model['weights'])))+len(ix)+len(pure)
    return {'status':'PASS','dimension':43,'points':total,'complete_families':len(C),
            'explicit_tiny_points':len(ix),'added_points':len(added),'deleted_points':len(deleted),
            'gain':len(added)-len(deleted),'tiny_pairs_checked':tiny_pairs,
            'tiny_family_angles_checked':len(ix)*len(C),'tiny_pure_angles_checked':len(ix)*len(pure),
            'complete_family_and_pure_checks':family_checks,'refined_body_bounds':refined,
            'arithmetic':'40D coefficient metric, checked int64 products, Python integers and isqrt cutoffs',
            'population_scope':'Inherited complete-family identities and P48n minimum/design inputs; no fresh large-shell enumeration',
            'points_sha256':hashlib.sha256(points.read_bytes()).hexdigest()}


def main():
    if not __debug__:raise SystemExit('Run without -O: assertions are part of the verifier.')
    parser=argparse.ArgumentParser(description=__doc__)
    here=Path(__file__).resolve().parent
    parser.add_argument('--root',type=Path,default=here.parent)
    parser.add_argument('--points',type=Path)
    parser.add_argument('--previous',type=Path,default=here/'previous_points.npz')
    parser.add_argument('--delta',type=Path,default=here/'delta.npz')
    args=parser.parse_args();points=args.points or args.root/'partial43/explicit_points_43.npz'
    print(json.dumps(check(args.root,points,args.previous,args.delta),indent=2))


if __name__=='__main__':main()
