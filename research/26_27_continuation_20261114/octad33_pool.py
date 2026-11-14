"""Complete Golay-octad midpoint family (3,3,1^6), with exact blockers."""
import argparse,itertools,json,time,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'common'))
from golay import golay24
from image_search import above

def main():
    ap=argparse.ArgumentParser();ap.add_argument('dim',type=int);ap.add_argument('--source',type=Path);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--cutoff',type=int,default=4);ap.add_argument('--family',choices=['33','75','aa','six','flat'],default='33');a=ap.parse_args();a.out.mkdir(exist_ok=False);t=time.monotonic();d=a.dim;root=(Path(__file__).resolve().parent/'data');src=a.source or root/f'PackingStar/25D-31D/new_bounds_20261111/D{d}_{201180 if d==26 else 209862}_coordinates.npz';z=np.load(src)
    if a.source:O=z['rational'].astype(np.int16);OB=z['radical'].astype(np.int16)
    elif d==26:O=2*np.r_[z['rational_points'],z['quadratic_A']];OB=2*np.r_[np.zeros_like(z['rational_points']),z['quadratic_B']]
    else:O=2*z['points'];OB=np.zeros_like(O)
    eq=(~np.any(O[:,24:],axis=1))&(~np.any(OB,axis=1));quad=np.flatnonzero(np.any(OB,axis=1));aux=np.flatnonzero(~eq&~np.any(OB,axis=1));known={row.tobytes() for row in O[~np.any(OB,axis=1)]};bulk={O[i,:24].tobytes():int(i) for i in np.flatnonzero(eq)};rng=np.random.default_rng(26002733);perm=rng.permutation(aux);F=O.astype(np.float32);G=golay24();octads=np.flatnonzero(G.sum(1)==8);signs=np.array([v for v in itertools.product([-1,1],repeat=8) if np.prod(v)==1],np.int16);pairs=np.array(list(itertools.combinations(range(8),2)));T=[]
    for i,j in itertools.combinations(range(d-24),2):
        for s,u in itertools.product([-4,4],repeat=2):
            tail=np.zeros(d-24,np.int16);tail[[i,j]]=[s,u];T.append(tail)
    T=np.array(T);nt=len(T);pool={};hist={};minimum=999;tested=0
    def body_groups():
        if a.family=='aa':
            ss=np.array(list(itertools.product([-4,4],repeat=3)),np.int16)
            for i in range(24):
                for j,k in itertools.combinations([v for v in range(24) if v!=i],2):
                    U=np.zeros((8,24),np.int16);V=U.copy();U[:,i]=ss[:,0];V[:,i]=ss[:,0];U[:,j]=ss[:,1];V[:,k]=ss[:,2];pu=np.array([bulk.get((2*u).tobytes(),-1) for u in U]);pv=np.array([bulk.get((2*v).tobytes(),-1) for v in V]);yield U+V,np.c_[pu,pv]
        elif a.family=='six':
            ss=np.array(list(itertools.product([-1,1],repeat=6)),np.int16)
            for oc in octads:
                pos=np.flatnonzero(G[oc]);bodies=[];parents=[]
                for i,j in pairs:
                    rest=np.delete(pos,[i,j]);U=np.zeros((64,24),np.int16);U[:,rest]=2*ss;U[:,pos[i]]=2;U[:,pos[j]]=2*np.prod(ss,axis=1);V=U.copy();V[:,pos[[i,j]]]*=-1;pu=np.array([bulk.get((2*u).tobytes(),-1) for u in U]);pv=np.array([bulk.get((2*v).tobytes(),-1) for v in V]);bodies.append(U+V);parents.append(np.c_[pu,pv])
                yield np.concatenate(bodies),np.concatenate(parents)
        elif a.family=='flat':
            S=1-2*G.astype(np.int16)
            for i,j in itertools.combinations(range(24),2):
                U=S.copy();V=S.copy();U[:,i]*=-3;V[:,j]*=-3;pu=np.array([bulk.get((2*u).tobytes(),-1) for u in U]);pv=np.array([bulk.get((2*v).tobytes(),-1) for v in V]);yield U+V,np.c_[pu,pv]
        elif a.family=='33':
            for oc in octads:
                pos=np.flatnonzero(G[oc]);V=np.zeros((128,24),np.int16);V[:,pos]=2*signs;pv=np.array([bulk.get((2*v).tobytes(),-1) for v in V]);bodies=[];parents=[]
                for i,j in pairs:
                    U=np.zeros_like(V);U[:,pos[i]]=4*signs[:,i];U[:,pos[j]]=4*signs[:,j];pu=np.array([bulk.get((2*u).tobytes(),-1) for u in U]);bodies.append(U+V);parents.append(np.c_[pu,pv])
                yield np.concatenate(bodies),np.concatenate(parents)
        else:
            S=1-2*G.astype(np.int16)
            for i in range(24):
                V=S.copy();V[:,i]*=-3;pv=np.array([bulk.get((2*v).tobytes(),-1) for v in V])
                for j in range(24):
                    if i==j:continue
                    U=np.zeros_like(V);U[:,i]=-4*S[:,i];U[:,j]=4*S[:,j];pu=np.array([bulk.get((2*u).tobytes(),-1) for u in U]);yield U+V,np.c_[pu,pv]
    group_count=759 if a.family in ('33','six') else 24*253 if a.family=='aa' else 276 if a.family=='flat' else 24*23;body_count=759*(128 if a.family=='33' else 64)*28 if a.family in ('33','six') else 24*253*8 if a.family=='aa' else 276*4096 if a.family=='flat' else 24*23*4096
    print(json.dumps(dict(family=a.family,groups=group_count,bodies=body_count,tails=nt,source=len(O),cutoff=a.cutoff)),flush=True)
    for number,(H,body_parents) in enumerate(body_groups()):
        parents=np.repeat(body_parents,nt,axis=0);C=np.c_[np.repeat(H,nt,axis=0),np.tile(T,(len(H),1))];counts=np.sum(parents>=0,axis=1).astype(np.int32);active=np.arange(len(C));tested+=len(C)
        for at in range(0,len(perm),128):
            if not len(active):break
            ids=perm[at:at+128];counts[active]+=(C[active].astype(np.float32)@F[ids].T>64).sum(1);active=active[counts[active]<=a.cutoff]
        if len(quad) and len(active):
            fail=above(C[active].astype(np.float32)@F[quad].T,C[active].astype(np.float32)@OB[quad].astype(np.float32).T);counts[active]+=fail.sum(1);active=active[counts[active]<=a.cutoff]
        if len(active):
            D=C[active].astype(np.float32)@F[aux].T;Q=above(C[active].astype(np.float32)@F[quad].T,C[active].astype(np.float32)@OB[quad].astype(np.float32).T) if len(quad) else np.zeros((len(active),0),bool)
            for loc,i in enumerate(active):
                ca=C[i];key=ca.tobytes()
                if key in known or key in pool:continue
                bs=np.unique(np.r_[parents[i][parents[i]>=0],aux[np.flatnonzero(D[loc]>64)],quad[np.flatnonzero(Q[loc])]]);assert len(bs)==counts[i]
                pool[key]=(ca,bs);hist[len(bs)]=hist.get(len(bs),0)+1
                if len(bs)<minimum:minimum=len(bs);print('BEST',json.dumps(dict(octad=number,cost=minimum,pool=len(pool),seconds=time.monotonic()-t)),flush=True)
                if not len(bs):
                    np.savez_compressed(a.out/'candidate.npz',rational=np.r_[O,ca[None,:]],radical=np.r_[OB,np.zeros((1,d),np.int16)],radicand=3,squared_norm=128,inner_product_limit=64,coordinate_denominator=2);(a.out/'positive.json').write_text(json.dumps(dict(count=len(O)+1,gain=1,added=ca.tolist()))+'\n');print('POSITIVE',len(O)+1,flush=True);return
        if number%20==0:print('progress',number,group_count,len(pool),minimum,round(time.monotonic()-t,1),flush=True)
    vals=list(pool.values());A=np.array([v[0] for v in vals],np.int16).reshape(-1,d);bs=[v[1] for v in vals];ptr=np.r_[0,np.cumsum(list(map(len,bs)))];idx=np.concatenate(bs) if bs else np.array([],int);np.savez_compressed(a.out/'pool.npz',candidates=A,blocker_indptr=ptr,blocker_indices=idx);report=dict(dimension=d,source=str(src),octads=len(octads),tested=tested,pool=len(A),minimum=minimum,costs=hist,seconds=time.monotonic()-t);(a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('FINISHED',json.dumps(report),flush=True)
if __name__=='__main__':main()
