"""Larger multiplicity spaces for real C2/C3 actions, no antipodal constraint."""
import json,time,resource
import numpy as np
from scipy.optimize import check_grad
import cyclic_asymmetric as ca
P=ca.P
original=ca.matrices

def matrices(d,q,freq):
    # freq supplies irreducible frequencies, and leading trivial coordinates
    # fill the remaining dimension. Frequency q/2 is -I on a two-plane.
    trivial=d-2*len(freq);assert trivial>=0
    out=[]
    for t in range(q):
        G=np.eye(d)
        for j,f in enumerate(freq):
            a=2*np.pi*f*t/q;c=np.cos(a);s=np.sin(a);k=trivial+2*j
            G[k:k+2,k:k+2]=[[c,-s],[s,c]]
        out.append(G)
    return np.array(out)
ca.matrices=matrices

def main():
    start=time.perf_counter();rng=np.random.default_rng(20261105)
    G=matrices(5,2,[1,1]);v=rng.normal(size=100)
    e=check_grad(lambda z:ca.objective(z,20,5,G,True)[0],lambda z:ca.objective(z,20,5,G,True)[1],v);assert e<1e-5
    jobs=[(5,2,20,[1,1],True,1000,55),(5,3,14,[1,1],False,1000,55),(6,2,37,[1,1],False,600,55),(6,3,25,[1,1,1],False,600,55),(7,2,63,[1,1,1],True,300,55),(7,3,42,[1,1,1],True,300,55),(9,3,102,[1,1,1,1],True,20,55),(10,3,171,[1,1,1,1,1],False,10,55)]
    result=[ca.run(*j,rng) for j in jobs]
    report={'random_seed':20261105,'gradient_finite_difference_error':float(e),'outcome':'Bounded nonantipodal continuous search; no nonexistence result; no certified candidate.','results':result,'baselines':{'5':40,'6':72,'7':126,'9':306,'10':510},'wall_seconds':time.perf_counter()-start,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (P/'flexible_results.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
