"""Exact real stabilizer-state census and two-orbit obstruction.

No record coordinates are inputs. All states are reconstructed from affine
binary supports and degree <=2 Boolean signs. Dense all-pairs work is avoided.
"""
from __future__ import annotations

if not __debug__:
    raise SystemExit("Exact checks require assertions; rerun without Python -O or -OO.")
from itertools import combinations, product
from pathlib import Path
import json, time, hashlib, resource, platform
import numpy as np

OUT=Path(__file__).resolve().parent

def gaussian(m,r):
    a=b=1
    for j in range(r):
        a*=2**(m-j)-1; b*=2**(r-j)-1
    return a//b

def rref_bases(m,r):
    for pivots in combinations(range(m),r):
        free=[(i,j) for i,p in enumerate(pivots) for j in range(p+1,m) if j not in pivots]
        for label in range(2**len(free)):
            rows=[1<<p for p in pivots]
            for k,(i,j) in enumerate(free):
                if label>>k&1: rows[i]|=1<<j
            yield rows

def span(rows):
    return [xor_sum(rows,j) for j in range(2**len(rows))]

def xor_sum(rows,label):
    z=0
    for i,row in enumerate(rows):
        if label>>i&1:z^=row
    return z

def signs(r):
    labels=np.arange(2**r,dtype=np.uint32)
    terms=[np.ones(2**r,dtype=np.uint32)]
    terms += [(labels>>i)&1 for i in range(r)]
    terms += [((labels>>i)&1)*((labels>>j)&1) for i in range(r) for j in range(i)]
    values=np.arange(2**len(terms),dtype=np.uint32)
    bits=np.zeros((len(values),2**r),dtype=np.uint8)
    for i,term in enumerate(terms):bits^=((values>>i)&1).astype(np.uint8)[:,None]*term.astype(np.uint8)
    return 1-2*bits.astype(np.int8)

def enumerate_states(m):
    """Rows int8; normalize each row by sqrt(support_size)."""
    n=2**m; layers=[]; ranks=[]; counts=[]
    for r in range(m+1):
        S=signs(r); blocks=[]; spaces=0
        for basis in rref_bases(m,r):
            spaces+=1; V=span(basis); seen=set()
            for a in range(n):
                if a in seen:continue
                support=[a^v for v in V];seen.update(support)
                X=np.zeros((len(S),n),dtype=np.int8);X[:,support]=S;blocks.append(X)
        X=np.concatenate(blocks); expected=gaussian(m,r)*2**(m-r)*2**(1+r+r*(r-1)//2)
        assert spaces==gaussian(m,r) and len(X)==expected
        assert np.all(np.count_nonzero(X,axis=1)==2**r)
        assert len(np.unique(X,axis=0))==len(X)
        layers.append(X);ranks.extend([r]*len(X));counts.append(expected)
    return np.concatenate(layers),np.asarray(ranks,dtype=np.uint8),counts

def exact_representative_checks(X,ranks,m):
    """Six orbit representatives under affine permutations/quadratic signs.

    Compare signed integer sums, then square inequalities only for positive sums.
    No floating point threshold is involved in any conflict count.
    """
    n=2**m; supports=2**ranks.astype(np.int64); report=[]
    for r in range(m+1):
        v=np.zeros(n,dtype=np.int16);v[:2**r]=1
        dots=X.astype(np.int16)@v
        positive=dots>0
        bad=positive & (4*dots.astype(np.int64)**2 > 2**r*supports)
        selfsame=np.all(X==v,axis=1);bad &= ~selfsame
        within=bad & (ranks%2==r%2)
        cross=bad & (ranks%2!=r%2)
        assert within.sum()==0
        degree=int(cross.sum());assert degree==2*(n-1)
        distribution={str(s):int((cross & (ranks==s)).sum()) for s in range(m+1)}
        lower=2*(2**r-1);upper=2**(r+1)*(2**(m-r)-1)
        assert degree==lower+upper
        report.append({'support_dimension':r,'within_sector_conflicts':0,'cross_conflicts':degree,'cross_by_support_dimension':distribution,'formula_lower':lower,'formula_upper':upper})
    return report

def main():
    started=time.perf_counter();X,ranks,counts=enumerate_states(5)
    checks=exact_representative_checks(X,ranks,5)
    even=X[ranks%2==0];even_ranks=ranks[ranks%2==0]
    assert len(even)==146880
    # Export exact integer rows and denominators. No numerical certificate claim.
    path=OUT/'bw32_reconstructed_control.npz'
    np.savez_compressed(path,integer_rows=even,support_dimension=even_ranks,
                        normalization='divide each row by 2**(support_dimension/2)')
    # Coordinate puncturing is not an isometry. Count only vectors whose removed
    # coordinates are zero; these are honest coordinate sections.
    rng=np.random.default_rng(20261104);sections=[]
    for n in range(25,32):
        best=0;best_removed=None
        for trial in range(32):
            removed=np.sort(rng.choice(32,32-n,replace=False))
            count=int(np.all(X[:,removed]==0,axis=1).sum())
            if count>best:best=count;best_removed=removed.tolist()
        sections.append({'dimension':n,'sampled_coordinate_sets':32,'best_raw_stabilizer_pool':best,'removed_coordinates':best_removed,'valid_kissing_subset_upper_bound':best})
    formula=[]
    for m in (4,5,6):
        layers=[gaussian(m,r)*2**(m-r)*2**(1+r+r*(r-1)//2) for r in range(m+1)]
        assert sum(layers[::2])==sum(layers[1::2])
        formula.append({'dimension':2**m,'support_layers':layers,'total_real_stabilizer_vectors':sum(layers),'one_sector':sum(layers[::2]),'exact_two_sector_independence_number':sum(layers[::2]),'bipartite_conflict_degree':2*(2**m-1)})
    report={'outcome':'No new record. Known Clifford/Barnes-Wall control reconstructed; exact obstruction closes unions/subsets of both real stabilizer sectors.',
      'random_seed':20261104,'control_status':'explicit known construction; algebraic family proof plus exact representative checks; not an independent exhaustive all-pairs certificate',
      'formula_results':formula,'enumerated_dimension':32,'enumerated_rows':len(X),'support_layer_counts':counts,'exact_representative_checks':checks,'coordinate_section_samples':sections,
      'control_file':path.name,'control_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
      'verified_repository_baselines':{'25':197584,'26':201136,'27':209807,'28':221988,'29':231697,'30':247972,'31':251188,'32':352800,'64':335496064},
      'sources':['https://arxiv.org/abs/math/0207186','https://arxiv.org/abs/math/0001038','https://arxiv.org/html/2404.17677v2','README.md verified-baseline table, retrieved 2026-11-04'],
      'scope':'Real stabilizer states only. Does not bound arbitrary Clifford orbits, arbitrary finite-group representations, deformed states, or kissing numbers.',
      'wall_seconds':time.perf_counter()-started,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'python':platform.python_version(),'numpy':np.__version__}
    (OUT/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
