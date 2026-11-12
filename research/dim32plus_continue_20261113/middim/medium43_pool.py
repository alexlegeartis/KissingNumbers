"""Exact finite exchange pool for the current 43D medium/tiny layer.

Every proposed point uses a supplied exact medium body and one of the rational
head templates in the current proof. Blockers are computed using integer
square-root cutoffs. This script does not change the incumbent package.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json,sys,math,time
import numpy as np
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).parent
W=ROOT/'verifications/improved/dim41-46-p48-heads/weyl43';P=W/'medium_split43'
sys.path.insert(0,str(W));import verify_families as vf

def exact_beta(A,B):
    assert int(np.max(np.abs(A)))*int(np.max(np.abs(B)))*A.shape[1]<2**53
    dots=A.astype(float)@B.astype(float).T
    assert np.array_equal(dots,np.rint(dots))
    ints=dots.astype(np.int64);assert not np.any(ints%144)
    return (ints//144).astype(np.int16)

def cutoff(za,na,zb,nb,D):
    return (math.isqrt(int(na)*int(nb))-2*sum(int(a)*int(b) for a,b in zip(za,zb)))//(D*D)

def build(family=59,maxblock=10,extra=False):
    start=time.time()
    g,residues,mu,tiny=vf.geometry()
    f=np.load(W/'configuration_43.npz');old=np.load(W/'explicit_points_43.npz')
    raw=np.load(W/'tiny_bodies.npz')['H'].astype(np.int64)
    delta=json.loads((P/'candidate.json').read_text());D=int(delta['head_denominator'])
    C=f['classes'];R=f['r3'];Z=f['head_numerators'].astype(object)*(D//int(f['head_denominator']))
    OZ=old['head_numerators'].astype(object)*(D//int(old['head_denominator']));I=old['body_indices'];pure=f['pure_directions'].astype(object)
    assert D%int(f['head_denominator'])==D%int(old['head_denominator'])==0
    H=[];J=[]
    for j in np.flatnonzero(R==12):
        hj=np.load(P/'medium_families'/f'family_{int(C[j])}.npz')['H'].astype(np.int64)
        H.extend(hj);J.extend((int(j),i) for i in range(495))
    H=np.array(H,dtype=np.int64);J=np.array(J,dtype=int)
    B=np.load(P/'medium_families'/f'family_{int(C[family])}.npz')['H'].astype(np.int64)
    current=[]
    dm=set(delta['delete']);dt=set(delta['delete_tiny'])
    for idx,(j,row) in enumerate(J):
        if j*495+row not in dm:current.append(dict(kind='medium',family=int(j),row=int(row),source_index=idx))
    nm=len(current)
    med_ids=np.array([p['source_index'] for p in current]);allB=[H[med_ids]];allZ=[Z[J[med_ids,0]]];allR=[np.full(nm,12)]
    tiny_ids=np.array([i for i in range(len(I)) if i not in dt])
    for i in tiny_ids:current.append(dict(kind='tiny',old_index=int(i),body_index=int(I[i])))
    allB.append(raw[I[tiny_ids]]);allZ.append(OZ[tiny_ids]);allR.append(np.full(len(tiny_ids),10))
    for j,p in enumerate(delta['points']+delta['new_tiny']):
        if 'family' in p:
            h=np.load(P/'medium_families'/f'family_{int(C[p["family"]])}.npz')['H'][p['row']];r=12
        else:h=raw[p['body_index']];r=10
        current.append(dict(kind='prior_insert',index=j,point=p));allB.append(np.array([h]));allZ.append(np.array([p['head_numerators']],dtype=object));allR.append(np.array([r]))
    allB=np.concatenate(allB);allZ=np.concatenate(allZ);allR=np.concatenate(allR)
    assert len(allB)==len(current)
    beta=exact_beta(B,allB);selfbeta=exact_beta(B,B)
    zheads=[]
    for p in delta['points']:
        z=tuple(p['head_numerators'])
        if z not in zheads:zheads.append(z)
    z=tuple(map(int,Z[family]));zheads.append(z)
    if extra:
        for zz in json.loads((HERE/f'medium43_extraheads_f{family}.json').read_text())['heads']:
            z=tuple(zz)
            if z not in zheads:zheads.append(z)
    hn=[12*D*D+sum(a*a for a in z) for z in zheads]
    oldn=[int(r)*D*D+sum(int(a)*int(a) for a in z) for r,z in zip(allR,allZ)]
    fn=[int(r)*D*D+sum(int(a)*int(a) for a in z) for r,z in zip(R,Z)]
    dc=((g['digits'][C]-g['digits'][int(C[family])])%3)@g['powers']
    universal=12+R-18+g['mu'][dc]//4
    vertices=[]
    for t,(z,n) in enumerate(zip(zheads,hn)):
        if any(int(universal[j])>cutoff(z,n,Z[j],fn[j],D) for j in np.flatnonzero(R!=12)):continue
        if any((v:=2*sum(int(a)*int(b) for a,b in zip(z,p)))>0 and v*v>n*sum(int(a)*int(a) for a in p) for p in pure):continue
        cut=np.array([cutoff(z,n,zz,nn,D) for zz,nn in zip(allZ,oldn)],dtype=np.int16)
        bad=beta>cut
        for row in np.flatnonzero(np.sum(bad,axis=1)<=maxblock):
            blockers=np.flatnonzero(bad[row]).tolist()
            # A source point already present is not a new proposal.
            same=np.flatnonzero(np.all(allB==B[row],axis=1))
            if any(tuple(map(int,allZ[i]))==tuple(z) for i in same):continue
            vertices.append(dict(row=int(row),head=t,blockers=blockers))
        print('head',t,'admitted',sum(v['head']==t for v in vertices),'min blockers',int(np.min(np.sum(bad,axis=1))),flush=True)
    paircut=np.array([[cutoff(z,n,zz,nn,D) for zz,nn in zip(zheads,hn)] for z,n in zip(zheads,hn)],dtype=np.int64)
    result=dict(status='Exact proposal pool; no new bound until exchange chosen and independently replayed',family=family,head_denominator=D,
      heads=[list(z) for z in zheads],normals=hn,vertices=vertices,current=current,maxblock=maxblock,seconds=time.time()-start)
    suffix='_extra' if extra else ''
    (HERE/f'medium43_pool_f{family}{suffix}.json').write_text(json.dumps(result))
    np.savez_compressed(HERE/f'medium43_pool_f{family}{suffix}.npz',selfbeta=selfbeta,paircut=paircut)
    print('POOL',len(vertices),'blocker histogram',np.unique([len(v['blockers']) for v in vertices],return_counts=True),'seconds',result['seconds'],flush=True)
    return result

if __name__=='__main__':build(int(sys.argv[1]) if len(sys.argv)>1 else 59,int(sys.argv[2]) if len(sys.argv)>2 else 10,'--extra' in sys.argv)
