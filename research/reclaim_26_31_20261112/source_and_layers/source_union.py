"""Rebuild or verify an injective blocker certificate for two-source unions.

If every point exclusive to the older repository code is matched injectively
to an incompatible PackingStar point, no subset of their union can exceed the
PackingStar cardinality. This conclusion does not rely on a solver tolerance
or on untested within-source angles. It is restricted to this literal union.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix, hstack, vstack
from scipy.sparse.csgraph import maximum_bipartite_matching
from source_paths import PACKINGSTAR_DEFAULT, packingstar_file, repository_file

def greater3(p, q, limit=64):
    z = limit-p
    return np.where(q >= 0, (z < 0) | (3*q*q > z*z),
                    (z < 0) & (z*z > 3*q*q))

def read(d, directory):
    ps_path = packingstar_file(directory, d)
    our_path = repository_file(d)
    with np.load(ps_path, allow_pickle=False) as f:
        if d == 26:
            den = int(f['rational_denominator'])
            assert den in (1, 2)
            x = f['rational_points'].astype(np.int16)*(2//den)
            ds = f['quadratic_denominators']
            assert np.all(ds == 1) and np.all(f['quadratic_radicands'] == 3)
            qa, qb = 2*f['quadratic_A'].astype(np.int16), 2*f['quadratic_B'].astype(np.int16)
        else:
            den = int(f['coordinate_denominator'])
            assert den in (1, 2)
            x = f['points'].astype(np.int16)*(2//den)
            qa = qb = np.empty((0, d), np.int16)
    with np.load(our_path, allow_pickle=False) as f:
        if d == 26:
            assert np.all(f['radical'] == 0)
            y = 2*f['rational'].astype(np.int16)
        else:
            den = int(f['coordinate_denominator'])
            assert den in (1, 2)
            y = f['points'].astype(np.int16)*(2//den)
    assert np.all((x.astype(np.int64)**2).sum(1) == 128)
    assert np.all((y.astype(np.int64)**2).sum(1) == 128)
    assert np.all((qa.astype(np.int64)**2+3*qb.astype(np.int64)**2).sum(1) == 128)
    assert np.all((qa.astype(np.int64)*qb.astype(np.int64)).sum(1) == 0)
    xd = {row.tobytes(): i for i, row in enumerate(x)}
    yd = {row.tobytes(): i for i, row in enumerate(y)}
    assert len(xd) == len(x) and len(yd) == len(y)
    common = xd.keys() & yd.keys()
    ai = np.array([i for key, i in xd.items() if key not in common], np.int32)
    bi = np.array([i for key, i in yd.items() if key not in common], np.int32)
    return x, qa, qb, y, ai, bi, len(common), ps_path, our_path

def check_injection(x, qa, qb, y, expected, our_ids, ps_ids):
    assert np.array_equal(expected, our_ids), 'Every exclusive old point must be covered'
    assert len(ps_ids) == len(our_ids) and len(np.unique(ps_ids)) == len(ps_ids)
    assert np.all((ps_ids >= 0) & (ps_ids < len(x)+len(qa)))
    rational = ps_ids < len(x)
    dots = np.sum(y[our_ids[rational]].astype(np.int64)*x[ps_ids[rational]].astype(np.int64), axis=1)
    assert np.all(dots > 64), 'A rational matched pair is not a conflict'
    qi = ps_ids[~rational]-len(x)
    rows = y[our_ids[~rational]].astype(np.int64)
    p = np.sum(rows*qa[qi].astype(np.int64), axis=1)
    q = np.sum(rows*qb[qi].astype(np.int64), axis=1)
    assert np.all(greater3(p, q)), 'An algebraic matched pair is not a conflict'

def run(d, directory, output, verify_only):
    t0 = time.time()
    x, qa, qb, y, ai, bi, common, ps_path, our_path = read(d, directory)
    certificate = output/f'matching{d}.npz'
    graph_edges = None
    if verify_only:
        with np.load(certificate, allow_pickle=False) as f:
            our_ids, ps_ids = f['our_indices'], f['packingstar_indices']
    else:
        a, b = x[ai], y[bi]
        assert np.max(np.abs(x), initial=0) <= 16 and np.max(np.abs(y), initial=0) <= 16
        assert d*16**2 < 2**24
        af = a.astype(np.float32)
        graph = vstack([csr_matrix(b[at:at+256].astype(np.float32)@af.T > 64)
                        for at in range(0, len(b), 256)], format='csr')
        if len(qa):
            p, q = b.astype(np.int32)@qa.astype(np.int32).T, b.astype(np.int32)@qb.astype(np.int32).T
            graph = hstack([graph, csr_matrix(greater3(p, q))], format='csr')
        graph_edges = int(graph.nnz)
        matched = maximum_bipartite_matching(graph, perm_type='column')
        assert np.all(matched >= 0), 'No saturating injection: this source union remains open'
        our_ids = bi
        ps_ids = np.empty(len(matched), np.int32)
        rational = matched < len(ai)
        ps_ids[rational] = ai[matched[rational]]
        ps_ids[~rational] = len(x)+matched[~rational]-len(ai)
        np.savez_compressed(certificate, our_indices=our_ids, packingstar_indices=ps_ids)
    check_injection(x, qa, qb, y, bi, our_ids, ps_ids)
    result = dict(dimension=d, packingstar_points=len(x)+len(qa), repository_points=len(y),
                  common_rational_points=common, repository_exclusive_points=len(bi),
                  packingstar_exclusive_points=len(ai)+len(qa), conflict_graph_edges=graph_edges,
                  matched_blockers=len(ps_ids), status='PASS exact injective blocker certificate',
                  maximum_feasible_cardinality_in_union=len(x)+len(qa),
                  scope='Literal two-source union only; no global or larger-family exclusion',
                  packingstar_sha256=hashlib.sha256(ps_path.read_bytes()).hexdigest(),
                  repository_sha256=hashlib.sha256(our_path.read_bytes()).hexdigest(),
                  matching_sha256=hashlib.sha256(certificate.read_bytes()).hexdigest(), seconds=time.time()-t0)
    if not verify_only:
        (output/f'union{d}.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--packingstar', type=Path, default=PACKINGSTAR_DEFAULT)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--dimensions', type=int, nargs='+', default=list(range(26, 32)))
    ap.add_argument('--verify-only', action='store_true')
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    for dimension in args.dimensions:
        run(dimension, args.packingstar, args.output, args.verify_only)
