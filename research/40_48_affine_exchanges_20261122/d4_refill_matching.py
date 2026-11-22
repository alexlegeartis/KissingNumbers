"""Check a matching from nine complete 44D refill families to distinct incumbent points."""
from pathlib import Path
import sys

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
CERT = HERE.parents[1]/"verifications/improved/dim44-d4-section"
sys.path.insert(0, str(CERT))
from moments import DET, K, verify_families


def main():
    # The parent minimum and the four calibrating shell counts are certified by
    # the main 44D verifier. Recheck the affine population identities used here.
    verify_families()
    with np.load(CERT/"ambient.npz") as data:
        gram, normals = data["gram"], data["normals"]
    with np.load(CERT/"families.npz") as data:
        selected = {(tuple(c),int(r)) for c,r in zip(data["labels"],data["radii"])}
        populations = data["affine_counts"] @ [1,0,0,0,134]
    adj = np.array(sp.Matrix(K).adjugate(), dtype=np.int64)
    inv6 = adj*6//DET
    assert np.array_equal(K @ inv6, 6*np.eye(4,dtype=int))
    labels = {tuple([0]*4)}
    for row in adj:
        labels = {tuple((np.array(v)+j*row) % DET) for v in labels for j in range(6)}
    labels = sorted(labels)
    assert len(labels) == DET
    index = {v:i for i,v in enumerate(labels)}
    opposite = [index[tuple(-np.array(v)%DET)] for v in labels]
    _, orbit = np.unique(np.minimum(opposite,np.arange(DET)),return_inverse=True)
    with np.load(HERE/"d4_refill_matching.npz",allow_pickle=False) as data:
        candidates, blockers = data["candidates"], data["blockers"]
    assert candidates.shape == blockers.shape == (9,4488,48)
    used = [set(),set()]
    for left,right in zip(candidates,blockers):
        pair = []
        for side,(raw,radius) in enumerate(zip((left,right),(13*DET//3,5*DET))):
            assert raw.dtype == np.int8
            points = raw.astype(np.int64)
            assert np.all(((points @ gram)*points).sum(1) == 6)
            inner = points @ gram @ normals.T
            assert np.all(inner == inner[0])
            v = inner[0]
            label = tuple((v @ adj) % DET)
            assert label not in used[side]
            used[side].add(label)
            assert 6*DET - v @ adj @ v == radius
            scaled = 6*points-(v @ inv6 @ normals)
            assert len(np.unique(scaled,axis=0)) == len(scaled)
            assert np.all(((scaled @ gram)*scaled).sum(1)*DET == 36*radius)
            if side == 0:
                assert (label,radius) not in selected
                # Check that v is a shortest normal representative, so the
                # norm-six population table counts this entire candidate family.
                pool = np.indices((13,)*4).reshape(4,-1).T-6
                representatives = pool[np.all((pool @ adj) % DET == label,axis=1)]
                assert min(((representatives @ adj)*representatives).sum(1)) == v @ adj @ v
                assert populations[orbit[index[label]]] == len(points)
            else:
                assert (label,radius) in selected
            pair.append(scaled)
        dots = ((pair[0] @ gram)*pair[1]).sum(1)
        assert np.all(dots > 0) and np.all(4*dots*dots > 156*180)
    print("Matched all 40,392 candidates to distinct incompatible incumbent points.")


if __name__ == "__main__":
    main()
