"""Fresh mixed-period real cyclic orbit optimization, then symmetry release.
No record coordinates are imported. All outputs are numerical point sets.
"""
from __future__ import annotations
from pathlib import Path
import argparse,hashlib,json,time,resource
import numpy as np
from scipy.optimize import minimize,check_grad
P=Path(__file__).resolve().parent

class BudgetExpired(Exception):pass

def powers(d,q,frequencies):
    fixed=d-2*len(frequencies)
    if fixed<0:raise ValueError('Too many real irreducibles')
    ans=[]
    for t in range(q):
        G=np.eye(d)
        for j,f in enumerate(frequencies):
            a=2*np.pi*f*t/q;c=np.cos(a);s=np.sin(a);i=fixed+2*j
            G[i:i+2,i:i+2]=[[c,-s],[s,c]]
        ans.append(G)
    return np.array(ans)

def full_screen(V):
    normerr=float(np.max(abs(np.sum(V*V,axis=1)-1)))
    maximum=-1.;violations=0;near_duplicates=0
    for start in range(0,len(V),128):
        A=V[start:start+128]@V.T
        for i,row in enumerate(A):
            part=row[start+i+1:]
            if len(part):
                maximum=max(maximum,float(part.max()))
                violations+=int(np.count_nonzero(part>.5+1e-10))
                near_duplicates+=int(np.count_nonzero(part>1-1e-10))
    return {'points':len(V),'dimension':V.shape[1],'maximum_cosine':maximum,'violating_pairs_at_1e-10':violations,'near_duplicate_pairs':near_duplicates,'maximum_norm_squared_error':normerr,'pairs_checked':len(V)*(len(V)-1)//2,'numerically_admissible':maximum<=.5+1e-10 and normerr<1e-10,'exact_certificate':False}

class Model:
    def __init__(self,d,q,freq,periods,rng,independent_frames=1,compressed=True):
        self.d=d;self.q=q;self.freq=freq;self.periods=np.array(periods,dtype=int);self.k=len(periods)
        self.compressed=compressed and independent_frames==1
        self.mask=np.ones((self.k,d));fixed=d-2*len(freq)
        for i,L in enumerate(periods):
            if q%L:raise ValueError('Period must divide group order')
            for j,f in enumerate(freq):
                if f*L%q:self.mask[i,fixed+2*j:fixed+2*j+2]=0
            if not self.mask[i].any():raise ValueError('Requested orbit has no fixed subspace')
        G=powers(d,q,freq);self.starts=np.r_[0,np.cumsum(periods)[:-1]]
        self.owner=np.repeat(np.arange(self.k),periods)
        self.phase=np.concatenate([np.arange(L) for L in periods])
        frames=[]
        for _ in range(independent_frames):
            Q,_=np.linalg.qr(rng.normal(size=(d,d)));frames.append(Q)
        # Conjugation is represented as an output rotation. Separate groups may
        # have independently rotated frames, so cross-group checks stay full.
        self.frames=np.array([frames[i%independent_frames] for i in range(self.k)])
        self.maps=np.einsum('pad,pbd->pab',self.frames[self.owner],G[self.phase])
        self.deadline=float('inf');self.best=None;self.evaluations=0
    def points(self,z):
        Z=z.reshape(self.k,self.d)*self.mask
        norms=np.linalg.norm(Z,axis=1)
        if np.min(norms)<1e-14:raise ValueError('Zero masked seed')
        X=Z/norms[:,None]
        V=np.einsum('pab,pb->pa',self.maps,X[self.owner])
        return X,V,norms
    def objective(self,z,power=2):
        if time.monotonic()>self.deadline:raise BudgetExpired
        X,V,norms=self.points(z);self.evaluations+=1
        if self.compressed:
            # Convert representatives into the common physical frame too.
            R=np.einsum('iab,ib->ia',self.frames,X)
            A=R@V.T;E=np.maximum(A-.5,0);E[np.arange(self.k),self.starts]=0
            W=self.periods[:,None]
            if power==2:C=W*E;loss=float(.5*np.sum(W*E*E))
            else:C=2e4*W*E**3;loss=float(5e3*np.sum(W*E**4))
            gradR=C@V
            gradX=np.einsum('iab,ia->ib',self.frames,gradR)
            gradV=C.T@R
        else:
            A=V@V.T;E=np.maximum(A-.5,0);np.fill_diagonal(E,0)
            if power==2:loss=float(.5*np.sum(E*E));gradV=2*E@V
            else:loss=float(5e3*np.sum(E**4));gradV=4e4*(E**3)@V
            gradX=np.zeros_like(X)
        back=np.einsum('pab,pa->pb',self.maps,gradV)
        np.add.at(gradX,self.owner,back)
        gradX=(gradX-X*np.sum(gradX*X,axis=1)[:,None])/norms[:,None]
        gradX*=self.mask
        if self.best is None or loss<self.best[0]:self.best=(loss,z.copy(),power)
        return loss,gradX.ravel()

def run_job(config,rng,seconds=40,seed_arrays=None):
    start=time.monotonic();d=config['d'];periods=config['periods']
    model=Model(d,config['q'],config['frequencies'],periods,rng,config.get('frames',1),config.get('compressed',True))
    model.deadline=start+seconds;trials=[];best=None
    # Different restarts and quartic polish alternate. Every selected best is
    # ranked by a full-pair maximum, rather than by incomparable objectives.
    for trial in range(config.get('restarts',6)):
        if time.monotonic()>model.deadline:break
        if seed_arrays is not None and trial==0:z=np.asarray(seed_arrays).reshape(model.k,d).copy()
        elif best is not None and trial%3:z=best[1]+rng.normal(size=(model.k,d))*(.03 if trial%3==1 else .15)
        else:z=rng.normal(size=(model.k,d))
        power=2 if trial%2==0 else 4;model.best=None
        try:
            result=minimize(model.objective,z.ravel(),args=(power,),jac=True,method='L-BFGS-B',options={'maxiter':1800,'maxls':25,'ftol':1e-12,'gtol':1e-8,'maxcor':12})
            z=result.x;status=int(result.status);nit=int(result.nit)
        except BudgetExpired:
            if model.best is None:break
            z=model.best[1];status='wall_budget';nit=None
        X,V,_=model.points(z);screen=full_screen(V)
        trials.append({'trial':trial,'power':power,'status':status,'iterations':nit,'screen':screen})
        if best is None or screen['maximum_cosine']<best[0]:best=(screen['maximum_cosine'],X.copy(),V.copy(),screen)
        print(json.dumps({'job':config['name'],'trial':trial,'power':power,'maximum_cosine':screen['maximum_cosine'],'violations':screen['violating_pairs_at_1e-10'],'seconds':time.monotonic()-start}),flush=True)
        if screen['numerically_admissible']:break
    if best is None:raise RuntimeError('No completed bounded trial')
    filename=config['name']+'.npz'
    np.savez_compressed(P/filename,seeds=best[1],vectors=best[2],periods=model.periods,frequencies=np.array(config['frequencies']),frames=model.frames)
    report={'config':config,'wall_budget_seconds':seconds,'seconds':time.monotonic()-start,'evaluations':model.evaluations,'trials':trials,'best_screen':best[3],'file':filename,'status':'numerical only; no exact certificate'}
    (P/(config['name']+'.json')).write_text(json.dumps(report,indent=2)+'\n')
    return report

def free_release(source,rng,seconds=60,name=None):
    V=np.load(P/source['file'])['vectors'];V=V+rng.normal(size=V.shape)*.004;V/=np.linalg.norm(V,axis=1)[:,None];d=V.shape[1];name=name or source['config']['name']+'_released'
    config={'name':name,'d':d,'q':1,'frequencies':[],'periods':[1]*len(V),'frames':1,'compressed':False,'restarts':4,'source_trial':source['file'],'mode':'all individual point vectors free; inherited symmetry released'}
    # Identity frame makes the input point cloud the actual starting seed cloud.
    # Model's random common rotation is harmless since all pairs are invariant.
    return run_job(config,rng,seconds,V)

def gradient_checks():
    rng=np.random.default_rng(9911);errors=[]
    for frames,compressed in [(1,True),(2,False)]:
        model=Model(6,5,[1,2],[5,5,1,1],rng,frames,compressed)
        z=rng.normal(size=(4,6)).ravel()
        for power in (2,4):
            e=check_grad(lambda x:model.objective(x,power)[0],lambda x:model.objective(x,power)[1],z)
            # Quartic gradients have larger absolute scale; relative norm is
            # also recorded and bounded against finite differences.
            scale=max(1.,np.linalg.norm(model.objective(z,power)[1]));relative=float(e/scale)
            if relative>2e-6:raise RuntimeError(f'Gradient mismatch {frames,power,relative}')
            errors.append({'frames':frames,'compressed':compressed,'power':power,'relative_gradient_error':relative})
    return errors

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--phase',choices=['group','release'],default='group');ap.add_argument('--seconds',type=float,default=40);ap.add_argument('--release-count',type=int,default=3);args=ap.parse_args()
    start=time.monotonic();rng=np.random.default_rng(2026110412);checks=gradient_checks()
    if args.phase=='group':
        configs=[
          {'name':'d12_c15_mixed842','d':12,'q':15,'frequencies':[1,2,3,5,6],'periods':[15]*55+[5]*3+[1]*2,'restarts':6},
          {'name':'d11_c21_mixed605','d':11,'q':21,'frequencies':[1,2,3,6,7],'periods':[21]*28+[7]*2+[3],'restarts':6},
          {'name':'d12_c7_fixed842','d':12,'q':7,'frequencies':[1,1,2,2,3],'periods':[7]*120+[1]*2,'restarts':6},
          {'name':'d11_c5_605','d':11,'q':5,'frequencies':[1,1,1,2,2],'periods':[5]*121,'restarts':6},
          {'name':'d12_c11_mixed842','d':12,'q':11,'frequencies':[1,2,3,4,5],'periods':[11]*76+[1]*6,'restarts':6},
          {'name':'d12_c13_845','d':12,'q':13,'frequencies':[1,2,3,4,5,6],'periods':[13]*65,'restarts':6},
          {'name':'d11_c11_605','d':11,'q':11,'frequencies':[1,2,3,4,5],'periods':[11]*55,'restarts':6},
          {'name':'d12_c5_rotated842','d':12,'q':5,'frequencies':[1,1,1,2,2],'periods':[5]*168+[1]*2,'frames':2,'compressed':False,'restarts':4},
        ]
        reports=[run_job(c,rng,args.seconds) for c in configs]
    else:
        original=json.loads((P/'group_results.json').read_text())['runs']
        ordered=sorted(original,key=lambda r:r['best_screen']['maximum_cosine'])
        selected=[]
        for d in (12,11):
            selected.append(next(r for r in ordered if r['best_screen']['dimension']==d))
        selected += [r for r in ordered if r not in selected and r['best_screen']['points']==842]
        selected=selected[:args.release_count]
        reports=[free_release(r,rng,args.seconds) for r in selected]
    result={'random_seed':2026110412,'phase':args.phase,'gradient_checks':checks,'runs':reports,'seconds':time.monotonic()-start,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'baselines':{'11':604,'12':841},'source_inputs':'No record coordinates read; random seeds only, or explicitly named output from this wave for release.','status':'Numerical search. Failed point sets are not kissing configurations. No new record without independent exact or interval certification.'}
    (P/(args.phase+'_results.json')).write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
