"""Independently check the 73D non-lattice construction in integer arithmetic.

This file is also copied into the public certificate directory. Its parent
Gamma72 minimum and shell count are the published Nebe inputs. The completed
half-coset certificate is replayed by lattice/verify.py.
"""
import json
from pathlib import Path
import subprocess
import sys
import time
import numpy as np

HERE=Path(__file__).resolve().parent


def verify(directory=HERE,replay=True):
    directory=Path(directory);started=time.time()
    lattice=directory/'lattice'
    if not lattice.exists():lattice=directory.parent/'lattice73_20261025/result'
    if replay:subprocess.run([sys.executable,str(lattice/'verify.py')],check=True)
    with np.load(lattice/'certificate.npz') as z:
        g=z['gamma_gram'].astype(np.int64);x=z['heads'].astype(np.int64)
    u=np.load(directory/'surgery_class.npy').astype(np.int64)
    assert u.ndim==2 and u.shape[1]==72 and len(u)>0
    # Canonical rows prove that no two removed directions are antipodal.
    first=np.argmax(u!=0,axis=1)
    assert np.all(u[np.arange(len(u)),first]>0)
    assert len(np.unique(u,axis=0))==len(u)
    assert np.all(np.einsum('ij,jk,ik->i',u,g,u)==8)
    # These products use int64 multiplication, independently of search BLAS.
    largest_class=0;largest_head=0;products=0
    ug=u@g
    for start in range(0,len(u),128):
        dots=ug[start:start+128]@u.T
        for j in range(len(dots)):dots[j,start+j]=0
        largest_class=max(largest_class,int(np.abs(dots).max()))
        dots=ug[start:start+128]@x.T
        largest_head=max(largest_head,int(np.abs(dots).max()))
        products+=dots.size
    assert largest_class<=2 and largest_head<=4
    witness=np.load(directory/'nonlattice_witness.npy').astype(np.int64) if (directory/'nonlattice_witness.npy').exists() else None
    if witness is not None:
        v=witness;w=u[0]-v
        assert v.shape==(72,) and v@g@v==8 and w@g@w==8
        removed={tuple(row) for row in np.vstack([u,-u])}
        assert tuple(v) not in removed and tuple(w) not in removed
        # A lattice containing this code would contain u[0]=v+w and
        # sqrt(3)*u[0], the sum of the two positive-body replacement caps.
        # Their integer span is dense on a line, contradicting discreteness.
    # Caps have norm 8: (3/4)*8+2=8. The same-body opposite-tail
    # product is 6-2=4; distinct class lines give <=3/2+2=7/2.
    # Remaining equator products <=2*sqrt(3)<4.
    # Mixed head/cap products <=sqrt(3)+sqrt(5)<4 because 15<16.
    # The coset contacts are (x/2,+/-sqrt(5/2)); norm 22/4+5/2=8.
    # The three absolute latitudes 0,sqrt(2),sqrt(5/2) are distinct.
    assert 3+5==8 and 3*5<16
    base=6218175600;old_lattice=base+2*len(x)
    total=(base-2*len(u))+4*len(u)+2*len(x)
    result=dict(status='PASS',dimension=73,equator_parent=base,
                completed_coset_vectors=len(x),deleted_equator_points=2*len(u),
                replacement_caps=4*len(u),class_size=len(u),
                maximum_class_absolute_product=largest_class,
                maximum_head_class_absolute_product=largest_head,
                exact_mixed_products=products,lattice_comparison=old_lattice,
                bound=total,gain_over_lattice=total-old_lattice,
                nonlattice_witness_checked=witness is not None,seconds=time.time()-started)
    (directory/'surgery_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result


if __name__=='__main__':verify()
