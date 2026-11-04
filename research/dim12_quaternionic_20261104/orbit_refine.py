"""Break sign-orbit symmetry in native-search near-misses, aiming for 842."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS']:
    os.environ[key]='1'
import argparse,hashlib,itertools,json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from threadpoolctl import threadpool_limits
from exact_cells import exact_frames
from search import ROOT,normalize,metrics,relax
from continuation import optimize


def make_code(ids):
    data=np.load(ROOT/'runs'/'orbit_pool'/'exact_representatives.npz')
    aa=data['A'][ids];bb=data['B'][ids]
    signs=np.array(list(itertools.product([-1,1],repeat=3))).repeat(4,axis=1)
    mixed=((aa[:,None,:]+np.sqrt(2)*bb[:,None,:])*signs[None,:,:]/4).reshape(-1,12)
    c,_=exact_frames();axial=np.concatenate([c.reshape(-1,4),-c.reshape(-1,4)])/2
    axes=[]
    for j in range(3):
        x=np.zeros((24,12));x[:,j*4:j*4+4]=axial;axes.append(x)
    return np.vstack([*axes,mixed])


def prune(x,count,rng):
    w=np.maximum(x@x.T-.5,0)**2;np.fill_diagonal(w,0)
    kept=np.ones(len(x),dtype=bool)
    for _ in range(count):
        degree=w@kept
        degree[~kept]=-1
        degree+=rng.uniform(0,1e-12,size=len(x))
        kept[np.argmax(degree)]=False
    for _ in range(20):
        degree=w@kept;delta=0;move=None
        for r in np.flatnonzero(~kept):
            gains=degree+w[r];gains[~kept]=-1
            v=int(np.argmax(gains));change=degree[r]-gains[v]
            if change<delta-1e-15:delta=change;move=(r,v)
        if move is None:break
        r,v=move;kept[r]=True;kept[v]=False
    return x[kept],np.flatnonzero(~kept).tolist()


def trial(spec):
    path,index,steps,rounds,out=spec;path=Path(path);out=Path(out);threadpool_limits(1)
    start=time.monotonic();seed=84884220261104+index;rng=np.random.default_rng(seed)
    lines=path.read_text().splitlines();ids=np.array([int(s) for s in lines[1].split()])
    assert len(ids)==97 and len(set(ids.tolist()))==97
    x=normalize(make_code(ids));x,removed=prune(x,6,rng);initial=metrics(x)
    best=x.copy();bm=initial;trace=[]
    for r in range(rounds):
        x,opt=relax(x,steps);m=metrics(x);trace.append(dict(round=r,**m,**opt))
        if m['energy']<bm['energy']:best=x.copy();bm=m;np.save(out/f'pool_{index:04d}.npy',best)
        if bm['max_cos']<.5-1e-10:break
        x=best.copy()
        g=np.maximum(x@x.T-.5,0);np.fill_diagonal(g,0)
        stress=np.sum(g*g,axis=1);top=np.argsort(stress)[-int(rng.choice([8,16,32])):]
        x[top]+=rng.normal(scale=float(rng.choice([.015,.04,.10])),size=(len(top),12));x=normalize(x)
    x=best.copy();mm=bm;minimax=x.copy()
    for p in [4,8,16,32,64]:
        x,opt=optimize(x,700,power=p);m=metrics(x)
        if m['max_cos']<mm['max_cos']:minimax=x.copy();mm=m
    np.save(out/f'pool_{index:04d}.npy',best);np.save(out/f'pool_{index:04d}_minimax.npy',minimax)
    r=dict(index=index,seed=seed,source=str(path.relative_to(ROOT)),
           source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
           removed_848_indices=removed,initial=initial,best=bm,minimax=mm,trace=trace,
           elapsed_seconds=time.monotonic()-start)
    (out/f'pool_{index:04d}.json').write_text(json.dumps(r,indent=2)+'\n');return r


def main():
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True)
    p.add_argument('--trials',type=int,default=256);p.add_argument('--workers',type=int,default=96)
    p.add_argument('--steps',type=int,default=4000);p.add_argument('--rounds',type=int,default=3)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    files=list(a.input.glob('candidate_*.txt'))
    files.sort(key=lambda p:(int(p.read_text().split()[0]),p.name))
    # Different 97-orbit sets are distinct starts even when their rounded costs agree.
    files=files[:a.trials]
    cfg=vars(a).copy();cfg['input']=str(a.input.resolve());cfg['out']=str(a.out.resolve());cfg['actual_trials']=len(files)
    cfg['sources_sha256']={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ['orbit_refine.py','search.py','continuation.py','kernels.py']}
    (a.out/'config.json').write_text(json.dumps(cfg,indent=2)+'\n');results=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(trial,(str(p.resolve()),i,a.steps,a.rounds,str(a.out.resolve()))) for i,p in enumerate(files)]):
            r=f.result();results.append(r)
            print(json.dumps({k:r[k] for k in ['index','initial','best','minimax','elapsed_seconds']}),flush=True)
            (a.out/'summary.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':main()
