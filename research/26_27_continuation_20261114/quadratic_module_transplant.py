"""Price complete modules under exact rotations over Q(sqrt(2))."""
import argparse,itertools,json,time,math
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import maximum_flow,breadth_first_order
from radical_radii import bsign,qsign

def main():
    ap=argparse.ArgumentParser();ap.add_argument('dim',type=int);ap.add_argument('module',type=Path);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--trials',type=int,default=256);ap.add_argument('--cutoff',type=int,default=12);ap.add_argument('--seed',type=int,default=269);a=ap.parse_args();a.out.mkdir(exist_ok=False);start=time.monotonic();r=(Path(__file__).resolve().parent/'data');d=a.dim;k=d-24;m=8+k;z=np.load(a.module) if a.module.suffix=='.npz' else None
    if z is not None and 'P' in z:A0=8*z['Q'];B0=4*z['P'];den0=3
    else:
        M=z['points'] if z is not None else np.loadtxt(a.module,dtype=np.int16);A0=np.zeros_like(M);B0=4*M;den0=1
    assert A0.shape[1]==m;z=np.load(r/f'PackingStar/25D-31D/new_bounds_20261111/D{d}_{201180 if d==26 else 209862}_coordinates.npz')
    if d==26:O=2*np.r_[z['rational_points'],z['quadratic_A']];OB=2*np.r_[np.zeros_like(z['rational_points']),z['quadratic_B']]
    else:O=2*z['points'];OB=np.zeros_like(O)
    I=np.r_[[0,2,4,5,6,10,11,23],np.arange(24,d)];P,inv,w=np.unique(np.c_[O[:,I],OB[:,I]],axis=0,return_inverse=True,return_counts=True);PA=P[:,:m].astype(np.float32);PB=P[:,m:].astype(np.float32);po=np.argsort(inv,kind='stable');pt=np.r_[0,np.cumsum(w)];rng=np.random.default_rng(a.seed);records=[];pool={};best=-99999
    def reflect(A,B,den):
        s=rng.choice([-1,1],8);s[-1]=np.prod(s[:-1]);da=A[:,:8].astype(np.int64)@s;db=B[:,:8].astype(np.int64)@s;sc=4//math.gcd(4,*map(int,da),*map(int,db));A=A*sc;B=B*sc;den*=sc
        for C in [A,B]:
            dot=C[:,:8].astype(np.int64)@s;assert not np.any(dot%4);C[:,:8]-=(dot//4)[:,None]*s
        return A,B,den
    for trial in range(a.trials):
        A=A0.astype(np.int64).copy();B=B0.astype(np.int64).copy();den=den0
        for _ in range(0 if trial<k+1 else (trial//(k+1))%3):A,B,den=reflect(A,B,den)
        core=rng.permutation(8);tail=8+rng.permutation(k);ct=trial%(k+1);pairs=list(zip(core[:ct],tail[:ct]));rest=core[ct:];rt=tail[ct:]
        count=len(rest)//2 if trial%4!=3 else int(rng.integers(len(rest)//2+1));pairs.extend((rest[2*j],rest[2*j+1]) for j in range(count));pairs.extend((rt[2*j],rt[2*j+1]) for j in range(len(rt)//2))
        for i,j in pairs:
            if np.any((A[:,i]+A[:,j])%2):A*=2;B*=2;den*=2
            ai=A[:,i].copy();aj=A[:,j].copy();bi=B[:,i].copy();bj=B[:,j].copy();A[:,i]=bi+bj;A[:,j]=bj-bi;B[:,i]=(ai+aj)//2;B[:,j]=(aj-ai)//2
        if trial>=k+1 and trial%3==0:A,B,den=reflect(A,B,den)
        div=math.gcd(den,*map(int,A.ravel()),*map(int,B.ravel()));A//=div;B//=div;den//=div
        if den>12:print('SKIP_DEN',trial,den,flush=True);continue
        assert np.all(np.sum(A*A+2*B*B,axis=1)==128*den*den) and not np.any(np.sum(A*B,axis=1));AF=A.astype(np.float32);BF=B.astype(np.float32);p=AF@PA.T;q=BF@PA.T
        bad=(bsign(p-64*den,q,AF@PB.T,BF@PB.T,2,3)>0) if d==26 else (qsign(p-64*den,q,2)>0);cost=bad@w
        for i in np.flatnonzero(cost<=a.cutoff):
            ca=np.zeros(d,np.int16);cb=np.zeros(d,np.int16);ca[I]=A[i];cb[I]=B[i];dv=math.gcd(den,*map(int,ca),*map(int,cb));ca//=dv;cb//=dv;dd=den//dv;rad=2 if cb.any() else 1;key=(rad,dd,np.r_[ca,cb].tobytes())
            if key not in pool:pool[key]=(ca,cb,dd,np.concatenate([po[pt[j]:pt[j+1]] for j in np.flatnonzero(bad[i])]))
        rr,cc=np.nonzero(bad);nm=len(A);npf=len(P);source=nm+npf;sink=source+1;row=np.r_[np.full(nm,source),rr,nm+np.arange(npf)];col=np.r_[np.arange(nm),nm+cc,np.full(npf,sink)];val=np.r_[np.ones(nm,np.int32),np.full(len(rr),nm+1,np.int32),w.astype(np.int32)];graph=coo_matrix((val,(row,col)),shape=(sink+1,sink+1)).tocsr();flow=maximum_flow(graph,source,sink);gain=nm-int(flow.flow_value);rec=dict(trial=trial,module=nm,gain=gain,denominator=den,min_cost=int(cost.min()),cheap_pool=len(pool),pairs=[list(map(int,v)) for v in pairs],seconds=time.monotonic()-start);records.append(rec)
        if gain>best or trial%10==0:best=max(best,gain);print(json.dumps(rec),flush=True)
        if gain>0:
            residual=graph-flow.flow;residual.data[residual.data<=0]=0;residual.eliminate_zeros();reach=breadth_first_order(residual,source,directed=True,return_predecessors=False);sel=reach[reach<nm];deleted=reach[(reach>=nm)&(reach<nm+npf)]-nm;keep=~np.isin(inv,deleted);CA=np.zeros((len(sel),d),np.int16);CB=np.zeros_like(CA);CA[:,I]=A[sel];CB[:,I]=B[sel];RA=np.r_[O[keep],CA];RB=np.r_[OB[keep],CB];rad=np.r_[np.full(keep.sum(),3),np.full(len(sel),2)];dens=np.r_[np.ones(keep.sum(),int),np.full(len(sel),den)];assert len(RA)==len(O)+gain;np.savez_compressed(a.out/'candidate.npz',rational=RA,radical=RB,radicands=rad,denominators=dens,squared_norm=128,inner_product_limit=64,coordinate_denominator=2);np.savez_compressed(a.out/'exchange.npz',added_A=CA,added_B=CB,added_denominator=den,removed_indices=np.flatnonzero(~keep));print('POSITIVE_REQUIRES_VERIFICATION',len(RA),flush=True);break
        if trial%10==0:(a.out/'progress.json').write_text(json.dumps(dict(best_gain=best,records=records),indent=2)+'\n')
    vals=list(pool.values());np.savez_compressed(a.out/'pool.npz',candidate_A=np.array([v[0] for v in vals],np.int16).reshape(-1,d),candidate_B=np.array([v[1] for v in vals],np.int16).reshape(-1,d),radicands=np.array([2 if v[1].any() else 1 for v in vals]),denominators=np.array([v[2] for v in vals]),blocker_indptr=np.r_[0,np.cumsum([len(v[3]) for v in vals])],blocker_indices=np.concatenate([v[3] for v in vals]) if vals else np.array([],int));report=dict(dimension=d,module=len(A0),best_gain=best,cheap_candidates=len(vals),records=records,seconds=time.monotonic()-start);(a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('FINISHED',json.dumps({k:v for k,v in report.items() if k!='records'}),flush=True)

if __name__=='__main__':main()
