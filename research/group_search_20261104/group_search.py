"""Finite unitary orbit discovery; no record coordinates as input.
Numerical search results are NOT certificates. Run verify.py for exact controls.
"""
import argparse, json, itertools
from pathlib import Path
import numpy as np
from scipy.optimize import milp, Bounds, LinearConstraint
from scipy.sparse import lil_matrix
OUT=Path(__file__).resolve().parent

def canonical(v):
    v=v/np.linalg.norm(v)
    j=np.flatnonzero(abs(v)>1e-7)[0]
    return v*np.conj(v[j])/abs(v[j])
def key(v):
    return tuple(np.round(np.r_[v.real,v.imag],7))
def orbit(seed, gens, limit=4000, projective=True):
    normalize = canonical if projective else lambda v: v/np.linalg.norm(v)
    seed=normalize(seed); vs=[seed]; seen={key(seed)}
    for v in vs:
        for g in gens:
            w=normalize(g@v); k=key(w)
            if k not in seen:
                seen.add(k); vs.append(w)
                if len(vs)>limit: return None
    return np.array(vs)

def weil(p=3, rank=2, parity=1):
    """Even/odd subrepresentations of finite Fourier/quadratic-phase action.
    Coordinate basis indexed by F_p^rank. All generators commute with x -> -x.
    """
    pts=list(itertools.product(range(p),repeat=rank)); idx={x:i for i,x in enumerate(pts)}
    N=len(pts); basis=[]; seen=set()
    for x in pts:
        if x in seen: continue
        y=tuple(-a%p for a in x); seen.update((x,y))
        if x==y:
            if parity==1:
                v=np.zeros(N); v[idx[x]]=1; basis.append(v)
        else:
            v=np.zeros(N); v[idx[x]]=1/np.sqrt(2); v[idx[y]]=parity/np.sqrt(2); basis.append(v)
    B=np.array(basis).T; z=np.exp(2j*np.pi/p); gs=[]
    for a in range(rank):
        F=np.zeros((N,N),complex)
        for x in pts:
            for b in range(p):
                y=list(x); y[a]=b
                F[idx[tuple(y)],idx[x]]=z**(b*x[a])/np.sqrt(p)
        gs.append(F * (-1j if p%4==3 else 1))
        gs.append(np.diag([z**(x[a]**2) for x in pts]))
    for a in range(rank):
        for b in range(a): gs.append(np.diag([z**(x[a]*x[b]) for x in pts]))
    ans=[B.T@g@B for g in gs]
    assert all(np.max(abs(g.conj().T@g-np.eye(len(basis))))<1e-12 for g in ans)
    return ans

def phase_lift(lines,m,offset=0):
    return np.concatenate([lines*np.exp(1j*(offset+2*np.pi*k/m)) for k in range(m)])
def maxcos(X,Y=None):
    same=Y is None
    if same: Y=X
    best=-1.
    for start in range(0,len(X),128):
        G=(X[start:start+128]@Y.conj().T).real
        if same:
            k=np.arange(len(G)); G[k,start+k]=-1
        best=max(best,float(G.max()))
    return best

def select_blocks(blocks, seconds=20, representatives=None):
    edges=[]
    for i,(name,X) in enumerate(blocks):
        for j in range(i):
            if maxcos(X if representatives is None else representatives[i],blocks[j][1])>0.5+1e-8: edges.append((i,j))
    A=lil_matrix((len(edges),len(blocks)))
    for k,(i,j) in enumerate(edges): A[k,i]=A[k,j]=1
    res=milp(-np.array([len(x) for _,x in blocks],float),integrality=np.ones(len(blocks)),
       bounds=Bounds(0,1),constraints=LinearConstraint(A.tocsr(),-np.inf,1),
       options={'time_limit':seconds})
    chosen=[] if res.x is None else np.flatnonzero(res.x>.5).tolist()
    return chosen, {'status':int(res.status),'message':res.message,'edges':len(edges),
        'best':sum(len(blocks[i][1]) for i in chosen),
        'upper_bound':None if getattr(res,'mip_dual_bound',None) is None else float(-res.mip_dual_bound)}

def discover(gens, trials=60, limit=500, seed=20261104):
    """Sample cyclic subgroup linear characters via simple eigenvalues.
    Joint subgroup characters additionally sampled using intersecting eigenspaces.
    Only 1D common eigenspaces accepted. No claim of subgroup exhaustiveness.
    """
    rng=np.random.default_rng(seed); d=gens[0].shape[0]; reps=[]; seen=set(); diagnostics=[]
    candidates=[('coordinate',np.eye(d,dtype=complex)[:,i]) for i in range(d)]
    words=[]
    for t in range(trials):
        U=np.eye(d,dtype=complex)
        word=rng.integers(len(gens),size=int(rng.integers(3,18))).tolist()
        for k in word: U=gens[k]@U
        words.append(U)
        vals,V=np.linalg.eig(U)
        for j,z in enumerate(vals):
            if sum(abs(vals-z)<1e-6)==1: candidates.append((f'cyclic:{word}:{j}',V[:,j]))
        # All eigenspaces of U intersect eigenspaces of another word, not only
        # simple eigenvalues. A common line supplies a character of <U,W>.
        if t and t%3==0:
            W=words[int(rng.integers(t))]; wvals,_=np.linalg.eig(W)
            for z in np.unique(np.round(vals,8)):
                for w in np.unique(np.round(wvals,8)):
                    _,s,vh=np.linalg.svd(np.vstack((U-z*np.eye(d),W-w*np.eye(d))))
                    if s[-1]<1e-6 and s[-2]>1e-5:
                        candidates.append((f'joint:{t}:{z}:{w}',vh[-1].conj()))
    for label,v in candidates:
        k=key(canonical(v))
        if k in seen: continue
        X=orbit(v,gens,limit)
        if X is None:
            diagnostics.append({'seed':label,'over_limit':True}); continue
        ks={key(x) for x in X}; seen.update(ks)
        if any(key(X[0]) in old for old,_ in reps): continue
        reps.append((ks,X)); diagnostics.append({'seed':label,'lines':len(X)})
    return [X for _,X in reps],diagnostics

def run(p,rank,parity,trials=60,limit=500,generators=None,tag=None):
    gs=weil(p,rank,parity) if generators is None else list(generators)
    d=gs[0].shape[0]
    assert all(g.shape==(d,d) and np.max(abs(g.conj().T@g-np.eye(d)))<1e-10 for g in gs)
    tag=tag or f"weil_p{p}_r{rank}_s{parity}"
    assert Path(tag).name==tag
    np.savez_compressed(OUT/f"{tag}_generators.npz",generators=np.array(gs))
    orbits,diag=discover(gs,trials,limit)
    blocks=[]; representatives=[]; inventory=[]
    # Canonical projective representatives must be phased individually; a
    # uniform regular m-gon on each line makes representative choices material.
    # For p-th cyclotomic action use 2p phases, plus square/hexagon alternatives.
    for i,lines in enumerate(orbits):
        X=orbit(lines[0],gs,limit=6*limit,projective=False)
        accepted=[]
        if X is None:
            inventory.append({"orbit":i,"lines":len(lines),"vector_orbit_over_limit":True})
            continue
        if maxcos(X[:1],X[1:])>.5+1e-8:
            inventory.append({"orbit":i,"lines":len(lines),"vectors":len(X),"accepted_phase_orders":[]})
            continue
        for m in sorted(set((1,2,3,4,6,2*p))):
            Y=phase_lift(X,m)
            _,ix=np.unique(np.round(np.c_[Y.real,Y.imag],7),axis=0,return_index=True)
            Y=Y[np.sort(ix)]
            R=X[:1]*np.exp(2j*np.pi*np.arange(m)/m)[:,None]
            # Each R row represents a full orbit under the same unitary action.
            # Exclude identical copies when comparing within the candidate block.
            dots=(R@Y.conj().T).real
            dots[dots>1-1e-8]=-1
            c=float(dots.max())
            if c<=.5+1e-8:
                accepted.append(m)
                for t in range(12):
                    blocks.append((f'o{i}:m{m}:phase{t}/24',Y*np.exp(2j*np.pi*t/24)))
                    representatives.append(R*np.exp(2j*np.pi*t/24))
        inventory.append({'orbit':i,'lines':len(lines),'vectors':len(X),'accepted_phase_orders':accepted})
    if blocks:
        chosen,res=select_blocks(blocks,representatives=representatives)
        X=np.concatenate([blocks[i][1] for i in chosen])
        res.update(selected=[blocks[i][0] for i in chosen],maxcos=maxcos(X),
           real_rank=int(np.linalg.matrix_rank(np.c_[X.real,X.imag],tol=1e-7)))
        assert res['maxcos']<=0.5+1e-8, 'compressed screening failed full check'
        assert len(np.unique(np.round(np.c_[X.real,X.imag],7),axis=0))==len(X)
        np.savez_compressed(OUT/f'{tag}.npz',vectors=X,**{f'orbit{i}':x for i,x in enumerate(orbits)})
    else: res={'best':0}
    result={'representation':tag,'random_seed':20261104,'p':p,'rank':rank,'parity':parity,'complex_dimension':d,'trials':trials,
        'line_limit':limit,'inventory':inventory,'seeds':diag,'result':res,
        'cross_checks':'one representative per transitive vector orbit against the other full orbit; final all-pairs numerical check',
        'certificate':'numerical screening only; phase grid and sampled subgroups not exhaustive'}
    (OUT/f'{tag}.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!='seeds'}),flush=True)
    return result
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--p',type=int,default=3);ap.add_argument('--rank',type=int,default=2);ap.add_argument('--parity',type=int,default=1);ap.add_argument('--trials',type=int,default=40);ap.add_argument('--limit',type=int,default=500)
    ap.add_argument('--generators',help='NPZ with unitary generators array (k,d,d); finite-group assumption is caller responsibility')
    ap.add_argument('--tag',help='output filename stem for custom generators')
    a=ap.parse_args()
    if a.generators and not a.tag: ap.error('--generators requires --tag')
    gs=None if a.generators is None else np.load(a.generators)['generators']
    run(a.p,a.rank,a.parity,a.trials,a.limit,gs,a.tag)
