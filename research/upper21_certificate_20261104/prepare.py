"""Verify automorphisms, compute moment orbits, and diagonalize a group-algebra element."""
import json,subprocess,time
import numpy as np
from scipy.linalg import hadamard,eigh
from model import *

start=time.monotonic()
if not (DATA/'symmetry.npz').exists():
    import pynauty
    graph=pynauty.Graph(1024,adjacency_dict={u:list(map(int,N[u])) for u in range(1024)},
                       vertex_coloring=[{0},set(range(1,1024))])
    result=pynauty.autgrp(graph)
    np.savez_compressed(DATA/'symmetry.npz',generators=np.array(result[0]),orbits=np.array(result[3]))
Z=np.load(DATA/'symmetry.npz');gens=Z['generators'].astype(np.int32)
for g in gens:
    assert g[0]==0 and len(np.unique(g))==1024
    assert np.array_equal(np.sort(g[N],axis=1),np.sort(N[g],axis=1))
    for b in (1<<np.arange(10)):
        assert np.array_equal(g[np.arange(1024)^b],g^g[b])
_,word_orbit=np.unique(Z['orbits'],return_inverse=True)
with (DATA/'generators.bin').open('wb') as f:
    np.array([1024,len(gens)],np.uint32).tofile(f);gens.astype(np.uint32).tofile(f)
subprocess.run([str(HERE/'orbits.exe'),str(DATA/'generators.bin'),str(DATA/'pair.bin'),str(DATA/'triple.bin')],check=True)
tables=[];representatives=[]
for name in ('pair','triple'):
    q=np.fromfile(DATA/(name+'.bin'),np.uint32);k=int(q[0]);table=q[1:1+1024**2].reshape(1024,1024)
    reps=q[1+1024**2:];assert len(reps)==k
    tables.append(table);representatives.append(reps)
pair,triple=tables;pr,tr=representatives
assert np.array_equal(triple,triple.T)
assert np.array_equal(triple,triple[np.arange(1024)[:,None],np.arange(1024)[:,None]^np.arange(1024)])
for g in gens:
    assert np.array_equal(triple,triple[g[:,None],g])
W=hadamard(1024).astype(np.int16)
F=np.array([W[:,word_orbit==i].sum(1) for i in range(word_orbit.max()+1)]).T
F=np.unique(F,axis=0)
# A symmetric group-algebra element commutes with every invariant moment matrix.
# Its eigenspaces therefore block diagonalize all such matrices simultaneously.
R=np.zeros((1024,1024));rng=np.random.default_rng(2121)
for g,c in zip(gens,rng.normal(size=len(gens))):
    R[np.arange(1024),g]+=c;R[g,np.arange(1024)]+=c
ev,Q=eigh(R,driver='evd')
cuts=np.r_[0,np.flatnonzero(np.diff(ev)>1e-7)+1,1024]
assert all(ev[b-1]-ev[a]<1e-8 for a,b in zip(cuts[:-1],cuts[1:]))
np.savez_compressed(DATA/'moments.npz',word_orbit=word_orbit,pair_orbit=pair,triple_orbit=triple,
                    pair_reps=pr,triple_reps=tr,fourier=F,Q=Q,eigenvalues=ev,cuts=cuts)
report=dict(automorphisms_verified=len(gens),word_orbits=int(word_orbit.max()+1),
            ordered_pair_orbits=len(pr),triple_orbits=len(tr),fourier_rows=len(F),
            psd_blocks=len(cuts)-1,maximum_block=int(np.diff(cuts).max()),
            block_sizes=dict(zip(*[list(map(int,a)) for a in np.unique(np.diff(cuts),return_counts=True)])),
            seconds=time.monotonic()-start)
(DATA/'preparation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
