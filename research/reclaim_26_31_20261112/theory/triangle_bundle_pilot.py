"""Price exact three-point midpoint fibers in fixed orthogonal tail planes."""
import argparse,collections,itertools,json,sys,time
from pathlib import Path
import numpy as np
from source_utils import ROOT,REPO,add_source_arguments,load_integer_source
sys.path.insert(0,str(REPO/'common'))
from golay import golay24

def sqsign(a,b,d=3):
    a,b=int(a),int(b)
    if b==0:return (a>0)-(a<0)
    if a==0:return (b>0)-(b<0)
    if (a>0)==(b>0):return (a>0)-(a<0)
    v=a*a-d*b*b
    return ((a>0)-(a<0))*((v>0)-(v<0))

def exceeds(p,q):
    # p + sqrt(3) q > 24 sqrt(2), exactly.
    return sqsign(p,q)>0 and sqsign(p*p+3*q*q-1152,2*p*q)>0

def thresholds(tails,i,j,phase):
    out=np.zeros(len(tails),dtype=np.int16)
    cache={}
    for ix,t in enumerate(tails):
        u,v=int(t[i])*phase,int(t[j])*phase
        if (u,v) not in cache:
            m=max(2*np.sqrt(3)*u,-np.sqrt(3)*u+3*v,-np.sqrt(3)*u-3*v)
            k=int(np.floor(24*np.sqrt(2)-m))+1
            def bad(p):return any(exceeds(p+c*v,b*u) for b,c in ((2,0),(-1,3),(-1,-3)))
            assert bad(k) and not bad(k-1)
            cache[u,v]=k
        out[ix]=cache[u,v]
    return out

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--save-best',type=int,default=1200)
    ap.add_argument('--output-dir',type=Path,required=True)
    add_source_arguments(ap);args=ap.parse_args()
    start=time.time()
    f,source_record,X,_=load_integer_source(31,args)
    X=X.astype(np.int16);out=args.output_dir;out.mkdir(parents=True,exist_ok=True)
    assert int(np.abs(X).max())*24<2**24
    tn=np.sum(X[:,24:].astype(np.int32)**2,axis=1)
    eqids=np.flatnonzero(tn==0); eq=X[eqids,:24]//2
    assert np.all(X[eqids,:24]%2==0)
    G=golay24();gw=[sum(int(v)<<j for j,v in enumerate(row)) for row in G]; gs=set(gw)
    # Full independent membership check for the retained equator.
    for z in eq:
        a=np.abs(z);nz=int(np.count_nonzero(z))
        if nz==2:assert np.count_nonzero(a==4)==2
        elif nz==8:
            assert np.count_nonzero(a==2)==8
            assert sum(int(v!=0)<<j for j,v in enumerate(z)) in gs
            assert np.count_nonzero(z<0)%2==0
        elif nz==24:
            j=int(np.argmax(a));assert a[j]==3 and np.count_nonzero(a==1)==23
            signs=z.copy();signs[j]*=-1
            assert sum(int(v<0)<<j for j,v in enumerate(signs)) in gs
        else:raise ValueError(('non-Leech equator',z))
    eqlookup={z.tobytes():int(i) for z,i in zip(eq,eqids)}
    full=X[(tn==32)&np.all(np.abs(X[:,:24])==2,axis=1),:24]//2
    currentheads=set(sum(int(v<0)<<j for j,v in enumerate(z)) for z in full)
    chosen={}; fixed=((0,1),(2,3),(4,5))
    for gi,g in enumerate(gw):
        for i,j in itertools.combinations(range(24),2):
            h=g^(1<<i)^(1<<j)
            if h in currentheads or (i,j) in fixed:chosen[h]=(gi,i,j)
    assert currentheads.issubset(chosen)
    codes=np.array(list(chosen),dtype=np.uint32)
    H=np.array([[1-2*((int(h)>>j)&1) for j in range(24)]for h in codes],dtype=np.int16)
    meta=np.array(list(chosen.values()),dtype=np.int16)
    parents=[]
    for gi,i,j in meta:
        ids=[]
        for p in (i,j):
            u=(1-2*G[gi].astype(np.int16));u[p]*=-3
            if u.tobytes() in eqlookup:ids.append(eqlookup[u.tobytes()])
        parents.append(ids)
    nc=np.array([len(p) for p in parents],dtype=np.int32)
    movingids=np.flatnonzero(tn!=0); moving=X[movingids]
    colors=[(2*r,2*r+1,p) for r in range(3) for p in (1,-1)]
    T=np.array([thresholds(moving[:,24:],*c) for c in colors],dtype=np.int16)
    counts=np.empty((len(H),6),dtype=np.int32)
    B=moving[:,:24].T.astype(np.float32)
    for a in range(0,len(H),64):
        dots=H[a:a+64].astype(np.float32)@B
        for c in range(6):counts[a:a+64,c]=np.count_nonzero(dots>=T[c],axis=1)+nc[a:a+64]
        if a%2048==0:print('priced',a,'/',len(H),'minimum',int(counts[:a+len(dots)].min()),'seconds',time.time()-start,flush=True)
    ranks=np.argsort(counts,axis=None)[:args.save_best]
    hi,co=np.unravel_index(ranks,counts.shape)
    ptr=[0]; blockers=[]
    for h,c in zip(hi,co):
        dots=H[h].astype(np.float32)@B
        ix=movingids[np.flatnonzero(dots>=T[c])].tolist()+parents[h]
        assert len(ix)==counts[h,c]
        blockers.extend(sorted(ix));ptr.append(len(blockers))
    np.savez_compressed(out/'triangle_bundle_pool.npz',heads=H,head_codes=codes,head_meta=meta,colors=np.array(colors),all_blocker_counts=counts,chosen_heads=hi,chosen_colors=co,blocker_ptr=np.array(ptr),blocker_ids=np.array(blockers),source=f.name)
    report={'source':str(f),'source_points':len(X),'heads':len(H),'current_heads':len(currentheads),'fixed_edges':fixed,'colors':colors,'bundle_count':counts.size,'minimum_blockers':int(counts.min()),'quantiles':{str(q):float(np.quantile(counts,q)) for q in (0,.001,.01,.1,.5,.9,1)},'at_most':{str(k):int(np.count_nonzero(counts<=k))for k in (0,1,2,3,4,8,16,32,64,128,256,512,1024)},'saved_bundles':len(hi),'saved_union_blockers':len(set(blockers)),'elapsed_seconds':time.time()-start,'scope':'Exact pricing of every one of the three points in each fixed-plane bundle against current incumbent; no bundle compatibility optimization yet.'}
    report['source']=f.name;report['source_sha256']=source_record['sha256']
    (out/'triangle_bundle_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()
