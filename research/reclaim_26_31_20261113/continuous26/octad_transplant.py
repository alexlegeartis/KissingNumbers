"""Move the complete valid 380-point module onto other Golay octads.

Exact profile multiplicities turn each two-configuration union into a
weighted bipartite vertex cover. A minimum cut gives the largest union code.
Only the source-defined module is moved; the complete outside is allowed to
lose individual projected profiles, charged at their actual multiplicities.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import maximum_flow,breadth_first_order
from rational_tail_pool import greater3


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--trials',type=int,default=128)
    ap.add_argument('--out',type=Path,required=True);ap.add_argument('--seed',type=int,default=20261113)
    ap.add_argument('--collect',action='store_true');ap.add_argument('--cutoff',type=int,default=4)
    a=ap.parse_args();root=Path(__file__).resolve().parents[3];src=root.parent/'PackingStar/25D-31D/new_bounds_20261111/D26_201180_coordinates.npz'
    z=np.load(src);A=np.vstack([z['rational_points'],z['quadratic_A']]).astype(np.int16)
    B=np.vstack([np.zeros_like(z['rational_points']),z['quadratic_B']]).astype(np.int16)
    normal=np.asarray([0,2,4,5,6,10,11,23]);J=np.r_[normal,[24,25]]
    mask=np.all(A[:,[i for i in range(24) if i not in normal]]==0,axis=1)&np.all(B[:,[i for i in range(24) if i not in normal]]==0,axis=1)
    YA=A[mask][:,J];YB=B[mask][:,J];assert len(YA)==380
    coremask=np.sum(((A[:,:24]!=0)|(B[:,:24]!=0)).astype(np.int64)*(1<<np.arange(24)),axis=1)
    octads=[int(s,16) for s in (root/'research/reclaim_26_31_20261112/theory/golay_octads.hex').read_text().split()]
    octads=sorted(octads,key=lambda o:int(((coremask&~o)==0).sum()))
    rng=np.random.default_rng(a.seed);rows=[];start=time.monotonic();pool={}
    oldset={np.r_[aa,bb].tobytes() for aa,bb in zip(A,B)} if a.collect else set()
    for trial in range(a.trials):
        o=octads[trial%len(octads)];positions=np.asarray([i for i in range(24) if (o>>i)&1]);perm=np.arange(8)
        if trial>=len(octads):perm=rng.permutation(8)
        J=np.r_[positions[perm],[24,25]]
        profiles,inverse,weights=np.unique(np.c_[A[:,J],B[:,J]],axis=0,return_inverse=True,return_counts=True)
        PA=profiles[:,:10].astype(np.int32);PB=profiles[:,10:].astype(np.int32)
        bad=greater3(YA.astype(np.int32)@PA.T+3*(YB.astype(np.int32)@PB.T),YA.astype(np.int32)@PB.T+YB.astype(np.int32)@PA.T,16)
        if a.collect:
            cheap=np.flatnonzero(bad@weights<=a.cutoff)
            order=np.argsort(inverse,kind='stable');offset=np.r_[0,np.cumsum(weights)]
            for k in cheap:
                ca=np.zeros(26,np.int16);cb=np.zeros(26,np.int16);ca[J]=YA[k];cb[J]=YB[k];key=np.r_[ca,cb].tobytes()
                if key in oldset or key in pool:continue
                pieces=[order[offset[j]:offset[j+1]] for j in np.flatnonzero(bad[k])]
                blockers=np.concatenate(pieces) if pieces else np.empty(0,np.int32)
                pool[key]=(ca,cb,blockers)
        rr,cc=np.where(bad);m=len(YA);n=len(profiles);s=m+n;t=s+1
        I=np.r_[np.full(m,s),rr,m+np.arange(n)];K=np.r_[np.arange(m),m+cc,np.full(n,t)]
        V=np.r_[np.ones(m,np.int32),np.full(len(rr),m+1,np.int32),weights.astype(np.int32)]
        graph=coo_matrix((V,(I,K)),shape=(t+1,t+1)).tocsr();flow=maximum_flow(graph,s,t)
        gain=m-int(flow.flow_value);row=dict(trial=trial,octad=hex(o),positions=J.tolist(),old_sector=int(((coremask&~o)==0).sum()),profile_count=n,conflicts=len(rr),gain=gain,seconds=time.monotonic()-start,pooled_new=len(pool))
        if gain>0:
            residual=graph-flow.flow;residual.data[residual.data<=0]=0;residual.eliminate_zeros()
            reachable=breadth_first_order(residual,s,directed=True,return_predecessors=False)
            sel=reachable[reachable<m];deleted=reachable[(reachable>=m)&(reachable<m+n)]-m
            keep=~np.isin(inverse,deleted);addA=np.zeros((len(sel),26),np.int16);addB=np.zeros_like(addA);addA[:,J]=YA[sel];addB[:,J]=YB[sel]
            FA=np.vstack([A[keep],addA]);FB=np.vstack([B[keep],addB]);assert len(FA)==len(A)+gain
            out=a.out.with_name(a.out.stem+f'_gain{gain}_trial{trial}.npz')
            np.savez_compressed(out,rational=FA,radical=FB,radicand=3,squared_norm=32,inner_product_limit=16)
            np.savez_compressed(out.with_name(out.stem+'_exchange.npz'),removed_indices=np.flatnonzero(~keep),added_A=addA,added_B=addB)
            row.update(removed=int((~keep).sum()),added=len(sel),candidate_file=str(out))
        rows.append(row)
        if trial%10==0 or gain>0:print(json.dumps(row),flush=True)
        a.out.write_text(json.dumps(dict(source=str(src),source_count=len(A),seed=a.seed,best_gain=max(r['gain'] for r in rows),scope='Exact finite union selection; a positive output still needs complete independent geometry replay',trials=rows),indent=2)+'\n')
    if pool:
        vals=list(pool.values());lengths=np.asarray([len(v[2]) for v in vals],np.int32)
        np.savez_compressed(a.out.with_suffix('.npz'),candidate_A=np.asarray([v[0] for v in vals]),candidate_B=np.asarray([v[1] for v in vals]),blocker_indices=np.concatenate([v[2] for v in vals]),blocker_indptr=np.r_[0,np.cumsum(lengths)])


if __name__=='__main__':main()
