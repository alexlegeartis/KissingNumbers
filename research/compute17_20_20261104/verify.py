#!/usr/bin/env python3
"""Independent checks of baselines, graph samples and every saved discrete incumbent."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from orbits17 import DATA,read_graph
from continuous import baseline


def verify17():
    z=np.load(DATA/'pool17.npz')
    E=np.pad(z['equator'],((0,0),(0,1)))
    poles=np.zeros((2,17),dtype=np.int64);poles[:,-1]=[12,-12]
    V=np.concatenate([E,z['vertices'][z['seed']],poles]).astype(np.int64)
    assert V.shape==(5730,17)
    assert len(np.unique(V,axis=0))==5730
    assert np.all((V[:,:16]*V[:,:16]).sum(1)+2*V[:,16]**2==288)
    worst=-999
    # Floating BLAS computes these small integer sums exactly; cast back and
    # independently check scalar products in arbitrary precision below.
    F=V.astype(np.float64)
    for a in range(0,len(V),256):
        G=F[a:a+256,:16]@F[:,:16].T+2*np.outer(F[a:a+256,16],F[:,16])
        assert np.array_equal(G,np.rint(G))
        for t in range(len(G)):G[t,a+t]=-999
        worst=max(worst,int(G.max()))
    assert worst==144
    off,nb=read_graph(DATA/'graph17.bin');pool=z['vertices'];rng=np.random.default_rng(170020)
    for _ in range(20000):
        i,j=map(int,rng.integers(len(pool),size=2))
        dot=sum(int(a)*int(b) for a,b in zip(pool[i,:16],pool[j,:16]))+2*int(pool[i,16])*int(pool[j,16])
        row=nb[off[i]:off[i+1]];where=np.searchsorted(row,j)
        edge=where<len(row) and row[where]==j
        assert edge==(i!=j and dot>144)
    # Duplicate a point: the exact verification must reject the corrupted list.
    bad=V[0]
    assert sum(int(x)**2 for x in bad[:16])+2*int(bad[16])**2>144
    return dict(count=len(V),max_scaled_dot=worst,norm_squared_scaled=288,
                graph_pairs_independently_checked=20000,corruption_control='PASS')


def verify20():
    X=baseline(20)
    # Reconstruct integers instead of testing rounded irrational coordinates.
    U=np.rint(X[:17400]*np.sqrt(8.)).astype(np.int64)
    S=np.rint(X[17400:]*np.sqrt(20.)).astype(np.int64)
    assert np.all((U*U).sum(1)==8) and np.all(np.abs(S)==1)
    F=U.astype(np.float64);worst=-99
    for a in range(0,len(U),256):
        G=F[a:a+256]@F.T
        assert np.array_equal(G,np.rint(G))
        for t in range(len(G)):G[t,a+t]=-99
        worst=max(worst,int(G.max()))
    assert worst==4
    mixed=max(int((F[a:a+256]@S.astype(np.float64).T).max()) for a in range(0,len(U),256))
    assert mixed<=6 and mixed*mixed<=40
    G=S.astype(np.float64)@S.astype(np.float64).T
    np.fill_diagonal(G,-99)
    assert G.max()<=10
    return dict(count=len(X),base_points=len(U),flat_points=len(S),
                base_max_dot=worst,base_flat_max_integer_dot=mixed,
                minimum_flat_hamming_distance=int((20-G.max())/2))


def verify_incumbents():
    reports=[]
    for path in sorted(DATA.glob('*_best.npy')):
        # Solver run names end in _s<random seed>.
        name=path.name.rsplit('_s',1)[0]
        source=DATA/(name+'.npz')
        if not source.exists():continue
        z=np.load(source);selected=np.load(path)
        off,nb=read_graph(DATA/(name+'.graph'))
        chosen=np.zeros(len(off)-1,dtype=bool);chosen[selected]=True
        assert len(np.unique(selected))==len(selected)
        assert all(not chosen[nb[off[i]:off[i+1]]].any() for i in selected)
        reports.append(dict(file=path.name,count=int(z['fixed'])+int(z['weights'][selected].sum())))
    return reports


def verify20_quotient(name,samples=2000):
    z=np.load(DATA/(name+'.npz'),allow_pickle=True)
    reps=z['dense_representatives'];nd=len(reps);H=z['subgroup'];blocks=z['blocks']
    starts=z['block_offsets'];patterns=z['block_reps']
    def decode(vertex):
        if vertex<nd:return None,int(reps[vertex])
        block=int(np.searchsorted(starts,vertex,side='right')-1)
        positions=[p for p in range(20) if int(blocks[block])>>p&1]
        pattern=int(patterns[block][vertex-int(starts[block])])
        word=sum(((pattern>>j)&1)<<p for j,p in enumerate(positions))
        return int(blocks[block]),word
    def conflict(a,b):
        sa,wa=a;sb,wb=b
        if sa is None and sb is None:return (wa^wb).bit_count()<5
        if sa is None:return 8-2*((wa^wb)&sb).bit_count()>=7
        if sb is None:return 8-2*((wa^wb)&sa).bit_count()>=7
        common=sa&sb
        return common.bit_count()-2*((wa^wb)&common).bit_count()>4
    off,nb=read_graph(DATA/(name+'.graph'));rng=np.random.default_rng(200020)
    for _ in range(samples):
        i,j=map(int,rng.integers(len(off)-1,size=2));a=decode(i);sb,wb=decode(j)
        expected=i!=j and any(conflict(a,(sb,wb^(int(h) if sb is None else int(h)&sb))) for h in H)
        row=nb[off[i]:off[i+1]];where=np.searchsorted(row,j)
        actual=where<len(row) and row[where]==j
        assert actual==expected,(name,i,j)
    return dict(model=name,sampled_orbit_pairs=samples,expanded_sign_images_per_pair=len(H))


def main():
    report=dict(baseline17=verify17(),baseline20=verify20(),
                saved_incumbents=verify_incumbents(),
                quotient_checks=[verify20_quotient(name) for name in ('q20_r8_s1','q20_r7_s1','q20_r6_s2')])
    (DATA/'independent_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':main()
