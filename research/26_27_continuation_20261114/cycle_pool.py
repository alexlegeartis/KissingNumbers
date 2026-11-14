"""Search whole contact cycles from multiple Leech automorphic images."""
import argparse,itertools,time,json,sys
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from image_search import above
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'common'))
from golay import golay24

def main():
    ap=argparse.ArgumentParser();ap.add_argument('dim',type=int);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--trials',type=int,default=32);ap.add_argument('--seed',type=int,default=79);ap.add_argument('--extra-cycles',type=Path);a=ap.parse_args();a.out.mkdir(exist_ok=False);start=time.monotonic();r=(Path(__file__).resolve().parent/'data');d=a.dim;z=np.load(r/f'PackingStar/25D-31D/new_bounds_20261111/D{d}_{201180 if d==26 else 209862}_coordinates.npz');O=2*np.r_[z['rational_points'],z['quadratic_A']] if d==26 else 2*z['points'];OB=2*np.r_[np.zeros_like(z['rational_points']),z['quadratic_B']] if d==26 else np.zeros_like(O);N=np.array([0,2,4,5,6,10,11,23]);J=np.setdiff1d(np.arange(24),N);norm=np.sum(O[:,:24].astype(np.int64)**2+3*OB[:,:24].astype(np.int64)**2,axis=1);normal=~np.any(O[:,J],axis=1)&~np.any(OB[:,J],axis=1);fix=(norm!=128)&((norm!=96)|normal);F=O[fix];FB=OB[fix];L=np.load(r/'leech.npy').astype(np.int16);positive=L[np.arange(len(L)),np.argmax(L!=0,axis=1)]>0;lk={v.tobytes():i for i,v in enumerate(L) if positive[i]};canon=np.array([lk[(v if positive[i] else -v).tobytes()] for i,v in enumerate(L)]);blocked=np.zeros(len(L),bool)
    for at in range(0,len(F),64):blocked|=above(2*L.astype(np.float32)@F[at:at+64,:24].astype(np.float32).T,2*L.astype(np.float32)@FB[at:at+64,:24].astype(np.float32).T).any(1)
    base=int((~blocked).sum()+len(F));assert base==(196828 if d==26 else 196934);T=[]
    for i,j in itertools.combinations(range(d-24),2):
        for t in [-1,1]:
            v=np.zeros(d-24,np.int16);v[[i,j]]=[4,4*t];T.append(v)
    T=np.array(T);nt=len(T);FC=F[:,:24].astype(np.float32);FT=F[:,24:].astype(np.float32);BC=FB[:,:24].astype(np.float32);cycles={};Gcode=golay24();rng=np.random.default_rng(a.seed)
    def seedcycles(d):
        g=np.load(r/f'parentgraph{d}/graph.npz');h=g['bodies'];es=np.flatnonzero(h[np.arange(len(h)),np.argmax(h!=0,axis=1)]>0);H=h[es];P=canon[g['parents'][es]];u=P[:,0];v=P[:,1];graph=coo_matrix((np.ones(2*len(P)),(np.r_[u,v],np.r_[v,u])),shape=(len(L),len(L))).tocsr();_,labels=connected_components(graph,directed=False);tail=g['tails'][es];palette=[]
        for i,j in itertools.combinations(range(d-24),2):
            for s in [-1,1]:
                t=np.zeros(d-24,np.int16);t[[i,j]]=[4,4*s];palette.append(t)
        look={v.tobytes():i for i,v in enumerate(palette)};color=np.array([look[v.tobytes()] for v in tail]);out=[]
        for c in np.unique(labels[P]):
            ee=np.flatnonzero(labels[u]==c);vv=np.unique(P[ee])
            if len(ee)!=4 or len(vv)!=4 or not np.any(L[vv][:,J],axis=1).all():continue
            assert len(np.unique(color[ee]//2))==1;out.append((H[ee],color[ee],L[vv]))
        return out
    seeds26=seedcycles(26);seeds27=seedcycles(27);assert len(seeds26)==544
    def add(seeds,perm,sign,original=False):
        admitted=0
        for H0,off0,P0 in seeds:
            H=H0[:,perm]*sign;PP=P0[:,perm]*sign;sg=np.where(H[np.arange(4),np.argmax(H!=0,axis=1)]>0,1,-1);H=H*sg[:,None];order=np.lexsort(H.T[::-1]);H=H[order];off=off0[order].copy();origin=int(off[0]);off=(off%2)^(origin%2);parent=[]
            for p in PP:
                if p[np.flatnonzero(p)[0]]<0:p=-p
                parent.append(lk[p.tobytes()])
            parent=np.array(sorted(parent));assert len(np.unique(parent))==4
            if blocked[parent].any():continue
            key=np.c_[H,off].astype(np.int16).tobytes()
            if key in cycles:
                if original:cycles[key][4]=origin
                continue
            allow=np.ones(nt,bool);pc=np.abs(H.astype(np.float32)@FC.T);qc=np.abs(H.astype(np.float32)@BC.T)
            for phase in range(nt):allow[phase]=not above(pc+np.abs(T[off^phase].astype(np.float32)@FT.T),qc).any()
            if not allow.any():continue
            if original:assert allow[origin]
            cycles[key]=[H,off,parent,allow,origin if original else -1];admitted+=1
        return admitted
    print('BASE',add(seeds26 if d==26 else seeds27,np.arange(24),np.ones(24,np.int16),True),'other',add(seeds27 if d==26 else seeds26,np.arange(24),np.ones(24,np.int16)),'total',len(cycles),flush=True);affine=[(aa,bb) for aa in sorted({pow(x,2,23) for x in range(1,23)}) for bb in range(23)];rng.shuffle(affine)
    for tr,(aa,bb) in enumerate(affine[:a.trials]):
        perm=np.r_[(aa*np.arange(23)+bb)%23,23];sign=1-2*Gcode[int(rng.integers(len(Gcode)))].astype(np.int16);ad=add(seeds26 if d==26 else seeds27,perm,sign);print('IMAGE',tr,aa,bb,'new',ad,'total',len(cycles),'seconds',round(time.monotonic()-start,1),flush=True)
    if a.extra_cycles:
        zz=np.load(a.extra_cycles);extra=list(zip(zz['H'],zz['offset'],zz['parent_vectors']));print('EXTRA',add(extra,np.arange(24),np.ones(24,np.int16)),'total',len(cycles),flush=True)
    vals=list(cycles.values());H=np.concatenate([v[0] for v in vals]).astype(np.int16);off=np.concatenate([v[1] for v in vals]);P=np.array([v[2] for v in vals]);allowed=np.array([v[3] for v in vals]);orig=np.array([v[4] for v in vals]);nc=len(vals);bodycomp=np.repeat(np.arange(nc),4);n=nt*nc;flags=np.zeros((nc,nc),np.uint8);HF=H.astype(np.float32)
    for at in range(0,len(H),64):
        D=np.abs(HF[at:at+64]@HF.T);ii,jj=np.nonzero(D>32);ii+=at;keep=bodycomp[ii]<bodycomp[jj];ii=ii[keep];jj=jj[keep];ci=bodycomp[ii];cj=bodycomp[jj];par=off[ii]^off[jj];dd=D[ii-at,jj];ff=(1<<par).astype(np.uint8)
        if d==27:ff|=np.where(dd>48,1<<(2+par),0).astype(np.uint8)
        ff|=np.where(dd>64,16,0).astype(np.uint8);np.bitwise_or.at(flags,(ci,cj),ff)
        if at%4096==0:print('GRAPH',at,len(H),'seconds',round(time.monotonic()-start,1),flush=True)
    owners={}
    for i,pp in enumerate(P):
        for p in pp:owners.setdefault(int(p),[]).append(i)
    for own in owners.values():
        if len(own)>1:
            ii,jj=np.array(list(itertools.combinations(own,2))).T;flags[ii,jj]|=16
    ci,cj=np.nonzero(flags);ff=flags[ci,cj];allowedflat=allowed.ravel();state=np.flatnonzero(allowedflat);where=np.full(n,-1,int);where[state]=np.arange(len(state));rr=[];cc=[]
    for i,j in itertools.combinations(range(nt),2):
        z=np.flatnonzero(allowed[:,i]&allowed[:,j]);rr.append((nt*z+i).astype(np.int32));cc.append((nt*z+j).astype(np.int32))
    for i,j in itertools.product(range(nt),repeat=2):
        delta=i^j;bad=(ff&16)!=0
        if i//2==j//2:bad|=(ff&(1<<delta))!=0
        if d==27:
            for p0 in range(2):
                if i//2!=j//2 or delta!=(p0^1):bad|=(ff&(1<<(2+p0)))!=0
        use=bad&allowed[ci,i]&allowed[cj,j];rr.append((nt*ci[use]+i).astype(np.int32));cc.append((nt*cj[use]+j).astype(np.int32))
    row=np.concatenate(rr);col=np.concatenate(cc);G=coo_matrix((np.ones(2*len(row),np.int8),(np.r_[where[row],where[col]],np.r_[where[col],where[row]])),shape=(len(state),len(state))).tocsr();G.data[:]=1;initial=where[nt*np.flatnonzero(orig>=0)+orig[orig>=0]];expected=544 if d==26 else 1616;assert len(initial)==expected and G[initial][:,initial].nnz==0
    with open(a.out/'graph.bin','wb') as f:
        np.array([len(state),len(initial)],'<u4').tofile(f);np.array([G.nnz],'<u8').tofile(f);initial.astype('<u4').tofile(f);G.indptr.astype('<u8').tofile(f);G.indices.astype('<u4').tofile(f)
    np.savez_compressed(a.out/'data.npz',H=H,offset=off,parents=P,allowed=allowed,state=state,initial=initial,fixed_A=F,fixed_B=FB,blocked=blocked,leech=L,base=base,dimension=d,tail_palette=T);report=dict(cycles=nc,states=len(state),initial=len(initial),directed_edges=int(G.nnz),seconds=time.monotonic()-start);(a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('FINISHED',json.dumps(report),flush=True)
if __name__=='__main__':main()
