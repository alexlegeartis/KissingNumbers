"""Translation and automorphism reduced theta-prime and lifted clique bounds.

All 21504 word-label candidates are retained. Averaging introduces no
periodicity restriction on an independent set.
"""
import argparse,json,subprocess,time
import numpy as np
import cvxpy as cp
from scipy.linalg import hadamard
from model import *


def prepare():
    path=DATA/'label_moments.npz'
    if path.exists():return np.load(path)
    gens=np.load(DATA/'symmetry.npz')['generators']
    with (DATA/'label_generators.bin').open('wb') as f:
        np.array([1024,len(gens)],np.uint32).tofile(f)
        gens.astype(np.uint32).tofile(f);P.astype(np.uint32).tofile(f)
    subprocess.run([str(HERE/'label_orbits.exe'),str(DATA/'label_generators.bin'),str(DATA/'label_orbits.bin')],check=True)
    b=np.fromfile(DATA/'label_orbits.bin',np.uint32);nv=int(b[0]);O=b[1:1+1024*441].reshape(1024,21,21);R=b[1+1024*441:]
    assert len(R)==nv and np.array_equal(O,O.transpose(0,2,1))
    edge=np.zeros(nv,dtype=bool)
    for d in range(1024):
        for a in range(21):
            for bb in range(21):
                bad=(d==0 and a!=bb) or (d in P and (a==bb or d==P[a] or d==P[bb]))
                if bad:edge[O[d,a,bb]]=True
    # Independently ensure each orbit has a uniform adjacency type.
    for d in range(1024):
        for a in range(21):
            for bb in range(21):
                bad=(d==0 and a!=bb) or (d in P and (a==bb or d==P[a] or d==P[bb]))
                assert edge[O[d,a,bb]]==bad
    W=hadamard(1024).astype(np.int16)
    dual=[]
    for g in gens:
        dual.append(np.array([sum(((int(c)&int(g[1<<i])).bit_count()%2)<<i for i in range(10)) for c in range(1024)]))
    covered=np.zeros(1024,bool);chars=[]
    for c in range(1024):
        if covered[c]:continue
        chars.append(c);queue=[c];covered[c]=True
        for a in queue:
            for g in dual:
                d=int(g[a])
                if not covered[d]:covered[d]=True;queue.append(d)
    blocks=[]
    for c in chars:
        blocks.append(np.array([np.bincount(O[:,a,b],weights=W[c],minlength=nv)
                                for a in range(21) for b in range(21)]).reshape(21,21,nv))
    row=np.bincount(O[:,0,:].ravel(),minlength=nv)
    cliques=[]
    for d in range(1024):
        base=np.bincount(O[d,0],minlength=nv)
        for p in range(21):
            r=base.copy();r[O[d^P[p],0,p]]+=1;cliques.append(r)
    cliques=np.unique(cliques,axis=0)
    np.savez_compressed(path,orbits=O,reps=R,blocks=blocks,characters=chars,
                        edges=edge,row=row,cliques=cliques)
    return np.load(path)


def four_cliques(z):
    path=DATA/'four_cliques.npy'
    if path.exists():return np.load(path)
    O=z['orbits'];nv=len(z['reps']);rows=np.zeros((1024*21*20,nv),dtype=np.int8);i=0
    for p in range(21):
        for q in range(21):
            if p==q:continue
            ds=np.arange(1024);ix=i+ds
            for words,label in ((ds,p),(ds,q),(ds^P[p],p),(ds^P[p],q)):
                np.add.at(rows,(ix,O[words,0,label]),1)
            i+=1024
    rows=np.unique(rows,axis=0)
    np.save(path,rows);return rows


def clique_products(z):
    """Moment inequalities for products of two maximal-clique slack indicators."""
    O=z['orbits'];nv=len(z['reps']);rows=[];sizes=[]
    for rep in z['reps']:
        d,p,q=map(int,(rep//441,(rep%441)//21,rep%21))
        base=np.bincount(O[d].ravel(),minlength=nv)
        left=np.bincount(O[d^P[p],p,:],minlength=nv)
        right=np.bincount(O[d^P[q],:,q],minlength=nv)
        cross=np.zeros(nv,dtype=np.int64);cross[O[d^P[p]^P[q],p,q]]=1
        for row,n in ((base,42),(base+left,43),(base+right,43),(base+left+right+cross,44)):
            rows.append(row);sizes.append(n)
    Q=np.unique(np.c_[rows,sizes],axis=0)
    return Q[:,:-1],Q[:,-1]


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--target',type=int);p.add_argument('--cliques',action='store_true');p.add_argument('--four',action='store_true');p.add_argument('--products',action='store_true');a=p.parse_args()
    start=time.monotonic();z=prepare();nv=len(z['reps']);x=cp.Variable(nv)
    con=[x>=0,x<=1,x[z['orbits'][0,0,0]]==1,x[np.flatnonzero(z['edges'])]==0]
    for block in z['blocks']:
        con.append(cp.reshape(block.reshape(441,nv)@x,(21,21),order='C')>>0)
    if a.cliques:
        con.append(z['cliques']@x<=1)
        if a.target:con.append(z['cliques']@x>=23-21504/a.target)
    if a.four:con.append(four_cliques(z)@x<=1)
    if a.products:
        assert a.target
        Q,s=clique_products(z)
        con.append(Q@x>=s-21504/a.target)
        con.append(Q@x<=np.floor(s/2))
    if a.target:
        con.append(z['row']@x==a.target);obj=cp.Minimize(0)
    else:obj=cp.Maximize(z['row']@x)
    prob=cp.Problem(obj,con)
    prob.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_feas=1e-9,tol_gap_rel=1e-9,max_iter=200)
    name='labels'+('_cliques' if a.cliques else '')+('_four' if a.four else '')+('_products' if a.products else '')+(f'_target{a.target}' if a.target else '_bound')
    report=dict(run=name,status=prob.status,objective=prob.value,target=a.target,
                variables=nv,psd_blocks=len(z['blocks']),block_size=21,clique_rows=len(z['cliques']),
                seconds=time.monotonic()-start,solver_seconds=prob.solver_stats.solve_time,numerical_only=True)
    if x.value is not None:
        np.save(DATA/(name+'_moments.npy'),x.value)
        report['maximum_clique_sum']=float((z['cliques']@x.value).max())
        report['minimum_psd_eigenvalue']=min(float(np.linalg.eigvalsh(np.einsum('abi,i->ab',b,x.value))[0]) for b in z['blocks'])
    (DATA/(name+'_result.json')).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
