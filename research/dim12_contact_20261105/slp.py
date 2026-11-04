"""Trust-region linearizations of the spherical minimax problem.

LP multipliers describe a local linear model. They are not a global upper bound
or a rigidity certificate. Every accepted step is checked on all actual pairs.
"""
from common import *
import argparse
import time
import warnings
from scipy.optimize import linprog, OptimizeWarning
from scipy.sparse import coo_matrix

def linear_model(x, radius, band=None, time_limit=180., threads=4):
    n, d = x.shape
    ii, jj = np.triu_indices(n, 1)
    dots = np.einsum('ij,ij->i', x[ii], x[jj])
    peak = float(dots.max())
    # An initial working set; missed constraints are harmless because all pairs
    # are checked after retraction. No global feasibility is inferred from it.
    if band is None:
        band = max(.002, 2 * np.sqrt(d) * radius)
    selected = dots > peak - band
    i, j, g = ii[selected], jj[selected], dots[selected]
    m = len(i)
    r = np.repeat(np.arange(m), 2*d+1)
    c = np.column_stack([i[:,None]*d+np.arange(d),
                         j[:,None]*d+np.arange(d), np.full(m,n*d)]).ravel()
    v = np.column_stack([x[j], x[i], -np.ones(m)]).ravel()
    aub = coo_matrix((v,(r,c)),shape=(m,n*d+1)).tocsr()
    aeq = coo_matrix((x.ravel(), (np.repeat(np.arange(n),d),np.arange(n*d))),
                     shape=(n,n*d+1)).tocsr()
    cost = np.zeros(n*d+1); cost[-1] = 1
    bounds = [(-radius,radius)]*(n*d) + [(None,None)]
    started = time.monotonic()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', OptimizeWarning)
        result = linprog(cost, A_ub=aub, b_ub=-g, A_eq=aeq, b_eq=np.zeros(n),
                         bounds=bounds, method='highs-ipm',
                         options={'time_limit':time_limit,'threads':threads,
                                  'dual_feasibility_tolerance':1e-8,
                                  'primal_feasibility_tolerance':1e-8,
                                  'ipm_optimality_tolerance':1e-10})
    info = dict(success=bool(result.success),status=int(result.status),message=result.message,
                radius=radius,band=band,constraints=m,initial_max=peak,
                elapsed_seconds=time.monotonic()-started)
    if not result.success:
        return None, None, info
    delta = result.x[:-1].reshape(n,d)
    lam = -result.ineqlin.marginals
    force = np.zeros_like(x)
    np.add.at(force,i,lam[:,None]*x[j]); np.add.at(force,j,lam[:,None]*x[i])
    force -= (force*x).sum(axis=1,keepdims=True)*x
    info.update(predicted_max=float(result.fun),predicted_decrease=peak-float(result.fun),
                weight_sum=float(lam.sum()),positive_weights=int((lam>1e-9).sum()),
                max_tangent_force=float(np.max(abs(force))),
                bounded_coordinates=int((abs(delta)>radius*(1-1e-6)).sum()),
                tangent_error=float(np.max(abs((delta*x).sum(axis=1)))))
    weights = dict(i=i,j=j,weight=lam,force=force,delta=delta)
    return delta, weights, info

def run(x, out, iterations=8, radius=.001, threads=4, time_limit=180.):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    x=normalize(x);trace=[]; best=metrics(x)
    np.save(out/'initial.npy',x)
    for step in range(iterations):
        delta,weights,info=linear_model(x,radius,threads=threads,time_limit=time_limit)
        info['iteration']=step
        if delta is None:
            trace.append(info);save_json(out/'trace.json',trace);print(json.dumps(info),flush=True);break
        np.savez(out/f'multipliers_{step:03d}.npz',**weights)
        candidates=[]
        for scale in [1.,.5,.25,.125,.0625]:
            y=normalize(x+scale*delta);m=metrics(y)
            candidates.append((m['max_cos'],scale,y,m))
        _,scale,y,m=min(candidates,key=lambda z:z[0])
        accepted=m['max_cos']<best['max_cos']-1e-12
        info.update(accepted=accepted,step_scale=scale,actual=m)
        if accepted:
            x=y;best=m
            radius=min(.02,radius*(1.3 if scale==1 else .7))
        else:
            radius*=.25
        trace.append(info);np.save(out/'best.npy',x);save_json(out/'trace.json',trace)
        print(json.dumps(info),flush=True)
        if best['max_cos']<=.5 or radius<1e-7:break
    save_json(out/'summary.json',dict(initial=metrics(np.load(out/'initial.npy')),best=best,trace=trace,
                                    source_sha256=sha(__file__)))
    return x,trace

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,default=OLD/'data'/'published842.npy')
    p.add_argument('--out',type=Path,default=ROOT/'runs'/'slp842');p.add_argument('--iterations',type=int,default=8)
    p.add_argument('--radius',type=float,default=.001);p.add_argument('--threads',type=int,default=4)
    p.add_argument('--time-limit',type=float,default=180.)
    a=p.parse_args();reserve();run(np.load(a.input),a.out,a.iterations,a.radius,a.threads,a.time_limit)
