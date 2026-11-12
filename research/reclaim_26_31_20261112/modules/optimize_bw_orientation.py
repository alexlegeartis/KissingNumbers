"""Bounded simulated annealing of a BW15 module against external obstructions."""
import argparse, hashlib, importlib.util, json, math, time
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
import bw15_transplant as bw

def mutate(A,kind,choice):
    Z=A.copy()
    if kind=='root':
        s=choice;v=Z[:,:8]@s
        if np.any(v%4):raise ValueError('E8 integrality')
        Z[:,:8]-=(v//4)[:,None]*s
    elif kind=='core_swap':
        i,j=choice;Z[:,[i,j]]=Z[:,[j,i]]
    elif kind=='tail_swap':
        i,j=choice;Z[:,[8+i,8+j]]=Z[:,[8+j,8+i]]
    else:Z[:,8+choice]*=-1
    return Z

def main():
    ap=argparse.ArgumentParser();ap.add_argument('source',type=Path);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--seed-trial',type=int,default=262);ap.add_argument('--iterations',type=int,default=20000)
    ap.add_argument('--seconds',type=float,default=150);ap.add_argument('--initial-module',type=Path)
    ap.add_argument('--rng-seed',type=int,default=20261113);ap.add_argument('--temperature',type=float,default=4)
    ap.add_argument('--restart-interval',type=int,default=1000)
    args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=False)
    start=time.monotonic();z=np.load(args.source);X=z['points'].astype(np.int16)
    if int(z['squared_norm'])==128:X//=2
    I=bw.NORMAL;J=[j for j in range(24) if j not in I]
    q=np.sum(X.astype(np.int64)[:,24:]**2,axis=1);h=np.sum(X[:,I].astype(np.int64)**2,axis=1)
    module=~np.any(X[:,J],axis=1);relevant=((q==8)&(h>8))|((q==16)&(h>0))
    ids=np.flatnonzero(relevant);external=relevant&~module
    proj=X[:,np.r_[I,np.arange(24,31)]];A0=proj[external]
    if args.initial_module:
        baseY=np.load(args.initial_module)['heads'].astype(np.int16)
        product=baseY.astype(np.int32)@baseY.astype(np.int32).T;np.fill_diagonal(product,-32)
        if baseY.shape!=(2016,15) or product.max()>16:raise ValueError('Initial module geometry')
    else:
        base=bw.bw15();r=np.random.default_rng(20261112)
        for t in range(args.seed_trial+1):baseY,_=bw.orientation(base,t,r)
    Yf=baseY.astype(np.float32);A=A0.copy();R=4*np.eye(15,dtype=np.int16)
    def cost(a):return int(np.any(Yf@a.astype(np.float32).T>16,axis=0).sum())
    score=cost(A);best=score;bestA=A.copy();bestR=R.copy();rng=np.random.default_rng(args.rng_seed);records=[]
    original={row.tobytes() for row in proj};best_gain=0
    def record_best(it):
        nonlocal best_gain
        numer=baseY.astype(np.int32)@R.astype(np.int32).T
        if np.any(numer%4):raise ValueError('Transformed module integrality')
        Y=(numer//4).astype(np.int16)
        # Exact transform and source projection identity.
        if not np.array_equal(A.astype(np.int64)*4,A0.astype(np.int64)@R.astype(np.int64)):raise ValueError('Transform identity')
        if not np.array_equal(R.astype(np.int32)@R.astype(np.int32).T,16*np.eye(15,dtype=np.int32)):raise ValueError('Isometry')
        C=np.asarray([row for row in Y if row.tobytes() not in original],np.int16)
        graph=csr_matrix(C.astype(np.float32)@proj[ids].astype(np.float32).T>16)
        matching=maximum_bipartite_matching(graph,perm_type='column');gain=len(C)-int(np.sum(matching>=0))
        rec={'iteration':it,'external_deletions':score,'whole_module_gain':58-score,'union_gain':gain,'new_points':len(C),'cross_conflicts':graph.nnz}
        records.append(rec);print(json.dumps(rec),flush=True)
        np.savez_compressed(args.out/'best_module.npz',heads=Y,transformation_numerator=R,transformation_denominator=4)
        if gain>0:
            reverse=np.full(len(ids),-1,np.int64)
            for i,j in enumerate(matching):
                if j>=0:reverse[j]=i
            left=set(map(int,np.flatnonzero(matching<0)));right=set();queue=list(left)
            for i in queue:
                for j0 in graph.indices[graph.indptr[i]:graph.indptr[i+1]]:
                    j=int(j0)
                    if matching[i]==j or j in right:continue
                    right.add(j)
                    if reverse[j]>=0 and int(reverse[j]) not in left:left.add(int(reverse[j]));queue.append(int(reverse[j]))
            removed=ids[sorted(right)];p=C[sorted(left)];added=np.zeros((len(p),31),np.int16);added[:,I]=p[:,:8];added[:,24:]=p[:,8:]
            keep=np.ones(len(X),bool);keep[removed]=False;F=np.r_[X[keep],added]
            if len(F)!=len(X)+gain or len(np.unique(F,axis=0))!=len(F):raise ValueError('Count/distinctness')
            for i,row in enumerate(added):
                dot=np.zeros(len(F),np.int32)
                for j in np.flatnonzero(row):dot+=F[:,j].astype(np.int32)*int(row[j])
                dot[len(X)-len(removed)+i]=-32
                if dot.max()>16:raise ValueError('Direct changed-pair check')
            np.savez_compressed(args.out/f'proposal31_{len(F)}.npz',points=2*F,squared_norm=128,inner_product_limit=64,coordinate_denominator=2)
            rec['exact_changed_pairs']='PASS';best_gain=max(best_gain,gain)
    record_best(0)
    completed=0
    for it in range(1,args.iterations+1):
        if time.monotonic()-start>=args.seconds:break
        u=rng.random()
        if u<.5:
            kind='root';choice=rng.choice([-1,1],8);choice[-1]=np.prod(choice[:-1])
        elif u<.65:kind='core_swap';choice=rng.choice(8,2,replace=False)
        elif u<.85:kind='tail_swap';choice=rng.choice(7,2,replace=False)
        else:kind='tail_sign';choice=int(rng.integers(7))
        candidate=mutate(A,kind,choice);value=cost(candidate);temperature=args.temperature*(1-(it%args.restart_interval)/args.restart_interval)+.1
        if value<=score or rng.random()<math.exp((score-value)/temperature):
            A=candidate;R=mutate(R,kind,choice);score=value
            if score<best:
                best=score;bestA=A.copy();bestR=R.copy();record_best(it)
        if it%args.restart_interval==0:A=bestA.copy();R=bestR.copy();score=best
        completed=it
    report={'source':str(args.source.resolve()),'seed_trial':args.seed_trial,'baseline':len(X),'iterations':completed,
        'initial_module':str(args.initial_module) if args.initial_module else None,
        'initial_module_sha256':hashlib.sha256(args.initial_module.read_bytes()).hexdigest() if args.initial_module else None,
        'rng_seed':args.rng_seed,'temperature':args.temperature,'restart_interval':args.restart_interval,
        'requested_iterations':args.iterations,'time_budget_seconds':args.seconds,'seconds':time.monotonic()-start,
        'best_external_deletions':best,'best_union_gain':best_gain,'records':records,
        'scope':'Bounded orientation search, no global optimum; whole-module and best-orientation finite-union counts only.'}
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':main()
