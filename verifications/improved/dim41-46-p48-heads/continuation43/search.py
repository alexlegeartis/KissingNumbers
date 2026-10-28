"""Explore continuous head positions after deleting a specified tiny-point set.

Numerical polytope vertices propose rational heads. Integer tests certify each
printed proposal against every retained point/family and pure head. Several
such individually feasible proposals still need mutual angle checks.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.spatial import HalfspaceIntersection

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent/'partial43'
sys.path.insert(0, str(SOURCE))
import verify_families as vf


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--body', type=int, required=True)
    parser.add_argument('--delete', nargs='*', type=int, default=[])
    args = parser.parse_args()
    g,residues,mu,_ = vf.geometry()
    X = np.load(SOURCE/'tiny_bodies.npz')['H'].astype(np.int64)
    assert 0 <= args.body < len(X)
    model = dict(np.load(SOURCE/'configuration_43.npz'))
    before = dict(np.load(HERE/'previous_points.npz'))
    assert all(0 <= i < len(before['body_indices']) for i in args.delete)
    body = np.load(SOURCE/'body_complement.npz')
    coefficients = vf.mm(vf.mm(X,body['b'].T),vf.inverse(body['q'],3))
    assert not np.any(coefficients%288)
    coefficients //= 288
    lookup = {row.tobytes():i for i,row in enumerate(residues)}
    labels = np.array([lookup[row.tobytes()] for row in coefficients%3])
    C,R = model['classes'],model['r3']
    dc = ((g['digits'][C]-g['digits'][labels[args.body]])%3)@g['powers']
    beta = 10+R-18+mu[dc]
    for j in np.flatnonzero(mu[dc]==8):
        pool = X[labels==dc[j]]
        if not np.any(vf.mm(X[args.body:args.body+1],pool.T)==(int(R[j])-20)*144):
            beta[j] -= 6
    old_indices = before['body_indices']
    raw = vf.mm(X[args.body:args.body+1],X[old_indices].T)[0]
    assert not np.any(raw%144)
    old_beta = raw//144
    Z = model['head_numerators'].astype(object)
    old_Z = before['head_numerators'].astype(object)
    d,old_d = int(model['head_denominator']),int(before['head_denominator'])
    pure = model['pure_directions'].astype(object)
    family_n = R.astype(object)*d*d+np.sum(Z*Z,axis=1)
    old_n = 10*old_d*old_d+np.sum(old_Z*old_Z,axis=1)
    pure_n = np.sum(pure*pure,axis=1)
    heads = np.r_[np.asarray(Z,dtype=float)/d,np.asarray(old_Z,dtype=float)/old_d]
    radii = np.r_[R,np.full(len(old_Z),10)]
    lengths = radii+np.sum(heads*heads,axis=1)
    planes = np.c_[np.r_[beta,old_beta]/np.sqrt(10*lengths),
                   2*heads/np.sqrt(lengths[:,None]),-np.ones(len(heads))]
    pure_float = np.asarray(pure,dtype=float)
    pure_planes = np.c_[np.zeros(len(pure)),2*pure_float/
                        np.linalg.norm(pure_float,axis=1)[:,None],-np.ones(len(pure))]
    planes = np.r_[planes,pure_planes,[[1,0,0,0,-1],[-1,0,0,0,0]]]
    keep = np.ones(len(planes),bool)
    keep[len(C)+np.array(args.delete,dtype=int)] = False
    vertices = HalfspaceIntersection(planes[keep],np.array([.01,0,0,0]),
                                     qhull_options='QJ1e-10').intersections
    D = 12000000000000
    proposals = []
    retained = np.setdiff1d(np.arange(len(old_Z)),args.delete)
    for v in vertices:
        if np.linalg.norm(v)<=1+1e-8 or v[0]<=1e-6:
            continue
        h = np.sqrt(10)*v[1:]/v[0]
        if abs(h).max()*D >= 2**62:
            continue
        z = np.rint(h*D).astype(np.int64).astype(object)
        N = 10*D*D+np.sum(z*z)
        value = beta.astype(object)*D*d+2*Z@z
        if not np.all((value<=0)|(value*value<=N*family_n)):
            continue
        value = old_beta[retained].astype(object)*D*old_d+2*old_Z[retained]@z
        if not np.all((value<=0)|(value*value<=N*old_n[retained])):
            continue
        value = 2*pure@z
        if np.all((value<=0)|(value*value<=N*pure_n)):
            proposals.append(tuple(map(int,z)))
    print(json.dumps(dict(body_index=args.body,deleted_indices=args.delete,
                         head_denominator=D,head_numerators=sorted(set(proposals)),
                         status='Each proposal individually passes exact tests; check mutual angles'),
                     indent=2))


if __name__ == '__main__':
    main()
