"""Mixed quadratic coefficients and individual denominators in exact exchanges."""
import argparse,json,math,time
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix,csr_matrix,vstack
from radical_radii import bsign,qsign

def main():
    ap=argparse.ArgumentParser();ap.add_argument('dim',type=int);ap.add_argument('pools',nargs='+',type=Path);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--cutoff',type=int,default=16);a=ap.parse_args();a.out.mkdir(exist_ok=False);start=time.monotonic();d=a.dim;r=(Path(__file__).resolve().parent/'data');z=np.load(r/f'PackingStar/25D-31D/new_bounds_20261111/D{d}_{201180 if d==26 else 209862}_coordinates.npz')
    if d==26:O=2*np.r_[z['rational_points'],z['quadratic_A']];OB=2*np.r_[np.zeros_like(z['rational_points']),z['quadratic_B']]
    else:O=2*z['points'];OB=np.zeros_like(O)
    known={(3 if b.any() else 1,1,np.r_[v,b].astype(np.int16).tobytes()) for v,b in zip(O,OB)};pool={}
    for path in a.pools:
        z=np.load(path);A=z['candidate_A' if 'candidate_A' in z else 'candidates'];B=z['candidate_B'] if 'candidate_B' in z else np.zeros_like(A);rr=z['radicands'] if 'radicands' in z else np.full(len(A),int(z['radicand']) if 'radicand' in z else 1 if not B.any() else 2);dd=z['denominators'] if 'denominators' in z else np.full(len(A),int(z['denominator']) if 'denominator' in z else 1);ptr=z['blocker_indptr'];idx=z['blocker_indices']
        for i in np.flatnonzero(np.diff(ptr)<=a.cutoff):
            ca=A[i].astype(np.int16);cb=B[i].astype(np.int16);rad=int(rr[i]) if cb.any() else 1;den=int(dd[i]);div=math.gcd(den,*map(int,ca),*map(int,cb));ca=ca//div;cb=cb//div;den//=div;key=(rad,den,np.r_[ca,cb].tobytes())
            if key in known:continue
            bs=np.sort(idx[ptr[i]:ptr[i+1]])
            if key in pool:assert np.array_equal(pool[key][-1],bs);continue
            pool[key]=(ca,cb,rad,den,bs)
    vals=list(pool.values());A=np.array([v[0] for v in vals],np.int16);B=np.array([v[1] for v in vals],np.int16);rad=np.array([v[2] for v in vals],np.int16);den=np.array([v[3] for v in vals],np.int16);bs=[v[4] for v in vals];old=np.unique(np.concatenate(bs));m=len(old);n=len(A);oldwhere={int(j):i for i,j in enumerate(old)};assert np.all(np.sum(A.astype(np.int64)**2+rad[:,None]*B.astype(np.int64)**2,axis=1)==128*den.astype(np.int64)**2);assert not np.any(np.sum(A.astype(np.int64)*B.astype(np.int64),axis=1));parts=[];AF=A.astype(np.float32);BF=B.astype(np.float32);print('LOADED',json.dumps(dict(candidates=n,old=m,fields={int(k):int(v) for k,v in zip(*np.unique(rad,return_counts=True))},denominators={int(k):int(v) for k,v in zip(*np.unique(den,return_counts=True))})),flush=True)
    for at in range(0,n,64):
        ids=np.arange(at,min(n,at+64));p=AF[ids]@AF.T;q=BF[ids]@AF.T;c=AF[ids]@BF.T;e=BF[ids]@BF.T;limit=64*den[ids,None].astype(np.int64)*den[None,:];same=rad[ids,None]==rad[None,:];bad=np.empty(p.shape,bool)
        if same.any():rr=np.broadcast_to(rad[ids,None],p.shape)[same];bad[same]=qsign(p[same]-limit[same]+rr*e[same],q[same]+c[same],rr)>0
        if (~same).any():rr=np.broadcast_to(rad[ids,None],p.shape)[~same];ss=np.broadcast_to(rad[None,:],p.shape)[~same];bad[~same]=bsign(p[~same]-limit[~same],q[~same],c[~same],e[~same],rr,ss)>0
        bad[np.arange(len(ids)),ids]=False;parts.append(csr_matrix(bad))
        if at%4096==0:print('GRAPH',at,n,round(time.monotonic()-start,1),flush=True)
    newG=vstack(parts,format='csr');assert (newG!=newG.T).nnz==0;ii,jj=newG.nonzero();up=ii<jj;ii=ii[up];jj=jj[up];br=np.array([oldwhere[int(j)] for x in bs for j in x],np.int32);bc=m+np.repeat(np.arange(n),[len(x) for x in bs]);row=np.r_[br,m+ii];col=np.r_[bc,m+jj];G=coo_matrix((np.ones(2*len(row),np.int8),(np.r_[row,col],np.r_[col,row])),shape=(m+n,m+n)).tocsr();keep=np.ones(len(O),bool);keep[old]=False;np.savez_compressed(a.out/'nodes.npz',A=np.r_[O[old],A],B=np.r_[OB[old],B],radicands=np.r_[np.where(np.any(OB[old],axis=1),3,1),rad],denominators=np.r_[np.ones(m,np.int16),den],common_A=O[keep],common_B=OB[keep],initial=np.arange(m),baseline_count=len(O),old_indices=old)
    with open(a.out/'graph.bin','wb') as f:
        np.array([m+n,m],'<u4').tofile(f);np.array([G.nnz],'<u8').tofile(f);np.arange(m,dtype='<u4').tofile(f);G.indptr.astype('<u8').tofile(f);G.indices.astype('<u4').tofile(f)
    report=dict(candidates=n,old=m,nodes=m+n,edges=int(G.nnz),seconds=time.monotonic()-start);(a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('FINISHED',json.dumps(report),flush=True)
if __name__=='__main__':main()
