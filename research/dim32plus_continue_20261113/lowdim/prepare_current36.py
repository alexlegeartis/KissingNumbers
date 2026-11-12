"""Exact phase admission for the current partially switched CT base.

Search preparation only. The accepted construction must be independently checked.
Each base triangle contains all three phases, so its phase-orbit admission
depends on two quotient labels, even when its large-vector switches are partial.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PACKAGE = ROOT / 'verifications/improved/dim36-37-coxeter-todd'
sys.path.insert(0, str(PACKAGE / 'phase'))
from exact_lattice import enumerate_short


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, default=HERE/'current36')
    ap.add_argument('--allow-lift-conflicts', action='store_true')
    ap.add_argument('--mode',choices=('equator','lifts'),default='equator')
    ap.add_argument('--equator-delta',type=Path)
    args = ap.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    start = time.time()
    source = PACKAGE/'completion/certificate.json'
    proposal = PACKAGE/'triangle_switch/proposal.json'
    c = json.loads(source.read_text())
    p = json.loads(proposal.read_text())
    G = np.array(c['gram_12'], dtype=np.int64)
    M = np.array(c['sublattice_matrix'], dtype=np.int64)
    phi = np.eye(12, dtype=np.int64)-M
    R0 = np.array(c['root_phase_zero_vectors'], dtype=np.int64)
    R = np.stack([R0, R0@phi, R0@phi@phi], axis=1)
    pow3 = 3**np.arange(6)
    first, second = np.array(c['label_first']), np.array(c['label_second'])
    def label(v): return ((v[:, first]+v[:, second])%3)@pow3
    short = list(enumerate_short(G.tolist()))
    V = np.array([v for v,n in short], dtype=np.int64)
    N = np.array([n for v,n in short], dtype=np.int64)
    L = label(V)
    W = np.full(729, 999, dtype=np.int64)
    W[0] = 0
    np.minimum.at(W, L, N)
    rootlabels = label(R0)
    assert np.array_equal(rootlabels, np.flatnonzero(W==4))
    assert len(V[N==8])==20412
    # The published triangle IDs use lexicographic large-vector order.
    Y = np.unique(V[N==8], axis=0)
    YL = label(Y)
    rootmap = {int(h):a for a,h in enumerate(rootlabels)}
    T = np.zeros((253,253,3), dtype=np.int8)
    for t in range(3):
        T[:252,:252,t] = (R0@G)@R[:,t].T
    T.tofile(out/'T.bin')
    RU = ((R.reshape(-1,12)@G)@Y.T).astype(np.int8).reshape(252,3,-1)
    # Maximizing a root phase against a fixed large vector is independent of
    # relative phase because every candidate includes its full three-point orbit.
    RY = RU.max(axis=1)
    del RU
    six = V[N==6]
    sixlabels = L[N==6]
    R6 = ((R0@G)@six.T).astype(np.int8)
    sixmax = {int(h):R6[:,sixlabels==h].max(axis=1)
              for h in np.flatnonzero(W==6)}
    # Full source six-cosets are phi-invariant. Check, rather than assume it.
    for h in sixmax:
        block=six[sixlabels==h]
        assert set(map(tuple,block))==set(map(tuple,block@phi))
    safe = {(i,j):np.ones((252,252), dtype=bool)
            for i,j in ((0,1),(0,2),(1,2))}
    def constrain(i,j,left,right):
        if i>j:
            i,j,left,right=j,i,right,left
        safe[i,j] &= left[:,None]+right[None,:] <= 6
    flipped=set(p['flipped_triangles'])
    owner=0
    for row in np.array(c['base_outer_cosets'],dtype=np.int64):
        norms=W[row]
        if sorted(norms)==[0,6,6]:
            i,j=np.flatnonzero(norms)
            constrain(int(i),int(j),sixmax[int(row[i])],sixmax[int(row[j])])
        elif sorted(norms)==[0,4,8]:
            i=int(np.flatnonzero(norms==4)[0])
            j=int(np.flatnonzero(norms==8)[0]);k=3-i-j
            a=rootmap[int(row[i])]
            yi=np.flatnonzero(YL==row[j]);assert len(yi)==81
            switch=np.array([owner*81+t in flipped for t in range(81)])
            small=T[:252,a].max(axis=1)
            if np.any(~switch):
                constrain(i,j,small,RY[:,yi[~switch]].max(axis=1))
            if np.any(switch):
                negative_small=(-T[:252,a].astype(np.int16)).max(axis=1)
                constrain(k,j,negative_small,RY[:,yi[switch]].max(axis=1))
            owner+=1
    assert owner==756
    # All root phases are compatible with the norm-twelve axial base.
    axial=V[N==4]@M
    assert int(((R0@G)@axial.T).max())<=6
    admissible=(safe[0,1][:,:,None]&safe[0,2][:,None,:]&safe[1,2][None,:,:])
    cells=np.argwhere(admissible).astype(np.int32)
    print('admitted cells',len(cells),'seconds',time.time()-start,flush=True)
    selected=set(map(tuple,c['phase_orbits']))
    selected.difference_update(map(tuple,p['remove_phase_orbits']))
    selected.update(map(tuple,p['add_phase_orbits']))
    assert len(selected)==32703
    if args.equator_delta:
        delta=json.loads(args.equator_delta.read_text())
        assert set(map(tuple,delta['remove_phase_orbits']))<=selected
        selected.difference_update(map(tuple,delta['remove_phase_orbits']))
        selected.update(map(tuple,delta['add_phase_orbits']))
    if args.mode=='lifts':
        pool=[]
        for k in range(3):
            i,j=[t for t in range(3)if t!=k]
            blocked=set()
            for orbit in selected:
                angles=(0,orbit[3],orbit[4])
                blocked.add((orbit[i],orbit[j],(angles[j]-angles[i])%3))
            for a,b in zip(*np.where(safe[i,j])):
                for relative in range(3):
                    if (a,b,relative)in blocked:continue
                    labels=[252]*3;labels[i]=int(a);labels[j]=int(b)
                    angles=[0]*3;angles[j]=relative
                    pool.append(labels+angles[1:])
        pool=np.array(pool,np.int32)
        seed=np.array(c['lift_phase_orbits'],np.int32)
        assert set(map(tuple,seed))<=set(map(tuple,pool))
        pool.tofile(out/'pool.bin');seed.tofile(out/'selected.bin')
        report={'mode':'lifts','current_base_switches':len(flipped),
                'equatorial_phase_orbits':len(selected),'pool_orbits':len(pool),
                'selected_orbits':len(seed),'maximum_body_product':2,
                'elapsed_seconds':time.time()-start}
        (out/'pool.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(report),flush=True);return
    phases=np.array([(u,v)for u in range(3)for v in range(3)],dtype=np.int32)
    pool=np.c_[np.repeat(cells,9,axis=0),np.tile(phases,(len(cells),1))]
    if not args.allow_lift_conflicts:
        forbidden=[np.zeros((253,253,3),dtype=bool)for _ in range(3)]
        for a,b,d,u,v in c['lift_phase_orbits']:
            labs=(a,b,d);angles=(0,u,v);k=labs.index(252)
            i,j=[t for t in range(3)if t!=k]
            forbidden[k][labs[i],labs[j],(angles[j]-angles[i])%3]=True
        bad=(forbidden[0][pool[:,1],pool[:,2],(pool[:,4]-pool[:,3])%3]
            |forbidden[1][pool[:,0],pool[:,2],pool[:,4]]
            |forbidden[2][pool[:,0],pool[:,1],pool[:,3]])
        pool=pool[~bad]
    seed=np.array(sorted(selected),dtype=np.int32)
    assert all(admissible[tuple(row[:3])]for row in seed)
    if not args.allow_lift_conflicts:
        assert not any(forbidden[0][b,d,(v-u)%3] or forbidden[1][a,d,v]
                       or forbidden[2][a,b,u] for a,b,d,u,v in seed)
    pool.tofile(out/'pool.bin');seed.tofile(out/'selected.bin')
    report={'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'proposal_sha256':hashlib.sha256(proposal.read_bytes()).hexdigest(),
            'base_triangles_switched':len(flipped),'base_points':403704,
            'admitted_cells':len(cells),'pool_orbits':len(pool),
            'selected_orbits':len(seed),'preserved_lifts':not args.allow_lift_conflicts,
            'safe_label_pairs':{str(k):int(v.sum())for k,v in safe.items()},
            'elapsed_seconds':time.time()-start}
    (out/'pool.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':main()
