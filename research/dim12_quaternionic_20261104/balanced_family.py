"""Three-parameter exact 840 family, and two simultaneous 25-point defects.

The M2,r weights are (1/2,sqrt(a_r),sqrt(3/4-a_r)), 1/4<=a_r<=1/2.
See BALANCED_FAMILY.md for the proof; the numerical search releases all symmetry.
"""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS']:
    os.environ[key]='1'
import argparse,hashlib,itertools,json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
from exact_cells import exact_frames,scalar_table
from search import ROOT,DATA,normalize,metrics,relax
from continuation import optimize


class Family:
    def __init__(self):
        c,d=exact_frames();c=c/2;d=d/np.sqrt(2)
        axial=np.concatenate([c.reshape(-1,4),-c.reshape(-1,4)])
        e1=np.zeros((24,12));e1[:,:4]=axial
        e2=np.zeros((24,12));e2[:,4:8]=axial
        e3=np.zeros((24,12));e3[:,8:]=axial
        pieces=[e1]
        for family in range(2):
            points=[]
            for r in range(3):
                delta=0 if family==0 else [2,1,3][r]
                for a,b in itertools.product(range(4),repeat=2):
                    cc=a^b^delta
                    for s in itertools.product([-1,1],repeat=3):
                        if family==0:q=[s[0]*d[r,a]/np.sqrt(2),s[1]*c[r,b]/2,s[2]*c[r,cc]/2]
                        else:q=[s[0]*c[r,a]/2,s[1]*d[r,b]/np.sqrt(2),s[2]*d[r,cc]/2]
                        points.append(np.concatenate(q))
            pieces.append(np.array(points))
        self.default_core=np.vstack(pieces);self.heads=np.vstack([e2,e3])
        self.group=np.repeat(np.arange(3),128)
        # Bring the published 841 displacement field into this explicit ordering.
        default=np.vstack([self.default_core,self.heads])
        source=np.load(DATA/'canonical840.npy');match=np.argmax(default@source.T,axis=1)
        assert len(set(match.tolist()))==840
        assert np.max(abs(default-source[match]))<1e-12
        record=np.load(DATA/'record841.npy')
        self.record=np.vstack([record[match],record[-1]])
        self.displacement=self.record[:792]-self.default_core
        self.swap=np.r_[np.arange(4),np.arange(8,12),np.arange(4,8)]
        balanced=self.core(np.ones(3)*.375)
        mapping=np.argmax(balanced[:,self.swap]@balanced.T,axis=1)
        assert len(set(mapping.tolist()))==792
        assert np.max(abs(balanced[:,self.swap]-balanced[mapping]))<1e-12
        self.swapped_displacement=np.empty_like(self.displacement)
        self.swapped_displacement[mapping]=self.displacement[:,self.swap]

    def core(self,a):
        x=self.default_core.copy();values=a[self.group]
        x[408:,4:8]*=np.sqrt(2*values)[:,None]
        x[408:,8:]*=np.sqrt(3-4*values)[:,None]
        return x


class HeadParameters:
    def __init__(self,family,n):self.family=family;self.n=n;self.calls=0
    def __call__(self,flat):
        self.calls+=1;n=self.n;raw=flat[:12*n].reshape(n,12);norms=np.linalg.norm(raw,axis=1,keepdims=True)
        h=raw/norms;a=flat[-3:];core=self.family.core(a)
        v=np.maximum(h@h.T-.5,0);np.fill_diagonal(v,0)
        cross=np.maximum(h@core.T-.5,0)
        value=.5*np.sum(v*v)+np.sum(cross*cross)
        gh=2*v@h+2*cross@core;gc=2*cross.T@h
        ga=np.zeros(3)
        for r in range(3):
            sl=slice(408+128*r,408+128*(r+1))
            ga[r]=(np.sum(gc[sl,4:8]*core[sl,4:8])/(2*a[r])
                   -np.sum(gc[sl,8:]*core[sl,8:])/(2*(.75-a[r])))
        gh-=np.sum(gh*h,axis=1,keepdims=True)*h
        return float(value),np.r_[(gh/norms).ravel(),ga]


def exact_support_check():
    c,d=exact_frames();table=scalar_table(c,d);checked=0
    for r,s in itertools.product(range(3),repeat=2):
        delta=[2,1,3][s]
        for a,b,aa,bb in itertools.product(range(4),repeat=4):
            cc=a^b;dd=aa^bb^delta
            nonzero=[bool(table[3+r,s,a,aa].any()),bool(table[r,3+s,b,bb].any()),
                     bool(table[r,3+s,cc,dd].any())]
            assert not all(nonzero);checked+=1
    return dict(exact=True,checked_m1_m2_label_pairs=checked,
                no_pair_has_three_nonzero_quaternionic_scalar_products=True,
                parameter_intervals=[[.25,.5]]*3,
                balanced_frame_eigenvalues=[66]*8+[78]*4)


def trial(spec):
    i,steps,rounds,out=spec;out=Path(out);threadpool_limits(1);seed=2537584220261104+i
    rng=np.random.default_rng(seed);start=time.monotonic();family=Family()
    if i%6==0:a=np.ones(3)*.375
    elif i%6==1:a=np.ones(3)*.35
    elif i%6==2:a=np.ones(3)*.40
    elif i%6==3:a=rng.uniform(.25,.5,size=3)
    elif i%6==4:a=rng.permutation([.3,.375,.45])
    else:a=np.clip(.375+rng.normal(scale=.02,size=3),.25,.5)
    if i%4==0:
        defect=np.vstack([family.record[816:840],family.record[-1]])
        heads=np.vstack([defect[:,family.swap],defect])
    elif i%4==2:
        defect=np.vstack([family.record[816:840],family.record[-1]])
        extra=np.zeros(12);extra[4:8]=rng.normal(size=4);extra/=np.linalg.norm(extra)
        heads=np.vstack([family.heads[:24],extra,defect])
    else:
        extra=np.zeros((2,12));extra[0,4:8]=rng.normal(size=4);extra[1,8:]=rng.normal(size=4)
        heads=np.vstack([family.heads,normalize(extra)])
    heads=normalize(heads+rng.normal(scale=[.005,.025,.06,.12][(i//4)%4],size=heads.shape))
    objective=HeadParameters(family,50)
    opt=minimize(objective,np.r_[heads.ravel(),a],jac=True,method='L-BFGS-B',
                 bounds=[(None,None)]*600+[(.25,.5)]*3,
                 options=dict(maxiter=1800,maxcor=12,ftol=1e-16,gtol=1e-10,maxls=30))
    heads=normalize(opt.x[:600].reshape(50,12));a=opt.x[-3:]
    core=family.core(a)
    structured=metrics(np.vstack([core,heads]))
    if i%3!=0:
        scale=[.25,.5,1.][(i//3)%3]
        core=normalize(core+scale*(family.displacement+family.swapped_displacement))
    x=np.vstack([core,heads]);best=x.copy();bm=metrics(x);trace=[]
    for r in range(rounds):
        x,info=relax(x,steps);m=metrics(x);trace.append(dict(round=r,**m,**info))
        if m['energy']<bm['energy']:best=x.copy();bm=m;np.save(out/f'balanced_{i:04d}.npy',best)
        if bm['max_cos']<.5-1e-10:break
        x=best.copy();chosen=rng.choice(np.arange(792,842),size=int(rng.choice([5,12,25,50])),replace=False)
        x[chosen]+=rng.normal(scale=float(rng.choice([.02,.05,.10])),size=(len(chosen),12));x=normalize(x)
    x=best.copy();mm=bm;minimax=x.copy()
    for power in [4,8,16,32]:
        x,info=optimize(x,600,power=power);m=metrics(x)
        if m['max_cos']<mm['max_cos']:minimax=x.copy();mm=m
    np.save(out/f'balanced_{i:04d}.npy',best);np.save(out/f'balanced_{i:04d}_minimax.npy',minimax)
    result=dict(index=i,seed=seed,parameters=a.tolist(),structured=structured,
                structured_calls=objective.calls,best=bm,minimax=mm,trace=trace,
                elapsed_seconds=time.monotonic()-start)
    (out/f'balanced_{i:04d}.json').write_text(json.dumps(result,indent=2)+'\n');return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--trials',type=int,default=512)
    p.add_argument('--workers',type=int,default=112);p.add_argument('--steps',type=int,default=4000)
    p.add_argument('--rounds',type=int,default=3);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    (ROOT/'balanced_family_exact_check.json').write_text(json.dumps(exact_support_check(),indent=2)+'\n')
    cfg=vars(a).copy();cfg['out']=str(a.out.resolve())
    cfg['sources_sha256']={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ['balanced_family.py','search.py','continuation.py','kernels.py','exact_cells.py']}
    (a.out/'config.json').write_text(json.dumps(cfg,indent=2)+'\n');results=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(trial,(i,a.steps,a.rounds,str(a.out.resolve()))) for i in range(a.trials)]):
            r=f.result();results.append(r)
            print(json.dumps({k:r[k] for k in ['index','parameters','structured','best','minimax','elapsed_seconds']}),flush=True)
            (a.out/'summary.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':main()
