"""Meaningful positive/negative verifier controls and independent numerical checks."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import itertools,json,tempfile
from pathlib import Path
import numpy as np
from kernels import NeighborList,pair_energy_gradient,power_energy_gradient
from continuation import PowerHinge
from orbit_search import OrbitHinge,group,generators
from exact_cells import find_clique
from search import Hinge,DATA,normalize
from verify import verify
from balanced_family import Family,HeadParameters,exact_support_check


def main():
    rng=np.random.default_rng(20261104);checks=[]
    # Exact equality and a violation much smaller than normal floating tolerances.
    with tempfile.TemporaryDirectory(prefix='dim12_checks_') as directory:
        path=Path(directory)/'control.npy'
        x=np.zeros((2,12));x[0,:2]=1;x[1,[0,2]]=1
        for name,delta,expected in [('equality',0.,True),('tiny_violation',-2**-40,False),('tiny_slack',2**-40,True)]:
            y=x.copy();y[1,2]+=delta;np.save(path,y)
            assert verify(path)['passed']==expected;checks.append(name)
        y=x.copy();y[1]=y[0];np.save(path,y);assert not verify(path)['passed'];checks.append('duplicate_rejected')
        for name,y in [('zero',np.zeros((2,12))),('nonfinite',np.full((2,12),np.nan)),('wrong_dimension',np.ones((2,13)))]:
            np.save(path,y)
            try:verify(path)
            except ValueError:checks.append(name+'_rejected')
            else:raise AssertionError(name)
    for dimension in [12,13,16]:
        x=normalize(rng.normal(size=(43,dimension)));neighbors=NeighborList(.5)
        for i in range(20):
            y=normalize(x+rng.normal(scale=.001 if i%4 else .1,size=x.shape))
            for power in [2,4,16]:
                v,g=neighbors.evaluate(y,power)
                vv,gg=(pair_energy_gradient(y,.5) if power==2 else power_energy_gradient(y,.5,power))
                assert abs(v-vv)<1e-11 and np.max(abs(g-gg))<1e-11
        checks.append(f'neighbor_bound_dimension_{dimension}')
    x=rng.normal(size=(19,13));v=rng.normal(size=x.shape);v/=np.linalg.norm(v);h=1e-6
    for power in [2,4,8,16,32]:
        fun=PowerHinge(x.shape,power=power,extra=1,penalty=.01);_,grad=fun(x.ravel())
        fd=(fun((x+h*v).ravel())[0]-fun((x-h*v).ravel())[0])/(2*h)
        assert abs(fd-grad@v.ravel())<1e-7
    checks.append('power_and_rank_penalty_gradients')
    rotations=group();z=generators(np.load(DATA/'canonical840.npy'),rotations);assert len(z)==280
    z=normalize(rng.normal(size=(15,12)));v=rng.normal(size=z.shape);v/=np.linalg.norm(v)
    f=OrbitHinge(z.shape,rotations);_,grad=f(z.ravel())
    fd=(f((z+h*v).ravel())[0]-f((z-h*v).ravel())[0])/(2*h)
    assert abs(fd-grad@v.ravel())<1e-7;checks.append('cyclic_group_and_gradient')
    # Compare the exhaustive clique search with a separate brute-force oracle.
    for _ in range(100):
        a=rng.uniform(size=(9,9))<.5;a=np.triu(a,1);a|=a.T
        bits=[sum(1<<int(j) for j in np.flatnonzero(row)) for row in a]
        actual=find_clique(bits,4)[0] is not None
        expected=any(all(a[i,j] for i,j in itertools.combinations(s,2)) for s in itertools.combinations(range(9),4))
        assert actual==expected
    checks.append('clique_search_vs_brute_force')
    family=Family();assert exact_support_check()['checked_m1_m2_label_pairs']==2304
    heads=normalize(family.heads+rng.normal(scale=.04,size=family.heads.shape))
    z=np.r_[heads.ravel(),np.ones(3)*.375];v=rng.normal(size=z.shape);v/=np.linalg.norm(v)
    f=HeadParameters(family,48);_,grad=f(z)
    fd=(f(z+h*v)[0]-f(z-h*v)[0])/(2*h)
    assert abs(fd-grad@v)<1e-7;checks.append('balanced_family_exact_support_and_gradient')
    import family_holes
    family_holes._family=family
    z=np.r_[normalize(rng.normal(size=(1,12)))[0],np.ones(3)*.375]
    v=rng.normal(size=z.shape);v/=np.linalg.norm(v);_,grad=family_holes.objective(z)
    fd=(family_holes.objective(z+h*v)[0]-family_holes.objective(z-h*v)[0])/(2*h)
    assert abs(fd-grad@v)<1e-7;checks.append('family_hole_gradient')
    result=dict(passed=True,checks=checks)
    (DATA.parent/'checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
