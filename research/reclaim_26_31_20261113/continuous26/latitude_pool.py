"""Joint exact rational-latitude candidate generation and finite MILP pilot."""
from pathlib import Path
import argparse,json,time,math
import numpy as np
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix
from rational_tail_pool import heads,greater3


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--seconds',type=int,default=120);ap.add_argument('--cutoff',type=int,default=8)
    ap.add_argument('--max-pool',type=int,default=6000)
    a=ap.parse_args();start=time.monotonic();z=np.load(Path(__file__).with_name('normal26.npz'))
    b=np.rint(z['B']).astype(np.int32);b=b[(b.astype(int)**2).sum(1)==16]
    old=np.rint(z['Y'][:364]).astype(np.int32)
    src=Path(__file__).resolve().parents[4]/'PackingStar/25D-31D/new_bounds_20261111/D26_201180_coordinates.npz'
    f=np.load(src);I=z['indices'];A=f['quadratic_A'][:,I].astype(int);B=f['quadratic_B'][:,I].astype(int)
    H=heads().astype(np.int64);D=65;tails=np.asarray([(x,y) for x in range(-D,D+1) for y in range(-D,D+1) if x*x+y*y==D*D],np.int64)
    specs=[(p,q) for q in [3,4,5,7,10,15,20] for p in [-1,1]]+[(0,1)]
    parts=[];dens=[];rows=[];seen={tuple(r) + (1,) for r in old}
    for p,q in specs:
        an=q*q-2*p*q-p*p;bn=q*q+2*p*q-p*p;d=q*q+p*p
        assert an*an+bn*bn==2*d*d
        ng=0;sel=[]
        for t in tails:
            v=np.zeros((len(H),10),np.int64);v[:,:8]=D*an*H;v[:,8:]=4*bn*t;den=d*D
            assert np.all(np.sum(v*v,axis=1)==32*den*den), 'Every generated point must have exact norm32'
            good=(v.astype(float)@b.astype(float).T<=16*den).all(1);v=v[good];ng+=len(v)
            bad=(v.astype(float)@old.T>16*den).sum(1)+greater3(v@A.T,v@B.T,16*den).sum(1)
            v=v[bad<=a.cutoff]
            for row in v:
                divisor=math.gcd(den,*map(int,row));rr=row//divisor;dd=den//divisor
                key=tuple(rr)+(dd,)
                if key not in seen:seen.add(key);sel.append((rr,dd))
        parts.extend(r[0] for r in sel);dens.extend(r[1] for r in sel)
        rows.append(dict(p=p,q=q,head_scale=[an,d],tail_scale=[bn,d],exterior_admissible=ng,retained=len(sel)))
        print(json.dumps(rows[-1]),flush=True)
    P=np.asarray(parts,np.int64).reshape((-1,10));pd=np.asarray(dens,np.int64)
    N=np.vstack([old,P]);nd=np.r_[np.ones(len(old),np.int64),pd];n=len(N)
    assert np.max(np.abs(N),initial=0)**2*10<2**53
    old_bad=np.c_[P.astype(float)@old.T>16*pd[:,None],greater3(P@A.T,P@B.T,16*pd[:,None])]
    groups=[[] for _ in range(380)];assignment=np.full(len(P),-1,np.int32)
    degree=old_bad.sum(1)
    for i in np.argsort(degree,kind='stable'):
        choices=np.flatnonzero(old_bad[i]);choices=sorted(choices,key=lambda j:len(groups[j]))
        for j in choices:
            members=groups[j]
            if not members or np.all(P[members].astype(float)@P[i].astype(float)>16*pd[members]*pd[i]):
                groups[j].append(int(i));assignment[i]=j;break
    missing=np.flatnonzero(assignment<0)
    cover_report=dict(cutoff=a.cutoff,latitude_specs=rows,new_rational_points=len(P),pool_points=len(P)+380,
                      scope='Exact integer conflict clique cover when every candidate is assigned',
                      covered=int(np.sum(assignment>=0)),uncovered=len(missing),baseline_sector=380,
                      seconds=time.monotonic()-start)
    np.savez_compressed(a.out.with_suffix('.npz'),rational=N,rational_denominators=nd,
                        quadratic_A=A,quadratic_B=B,assignment=assignment,indices=I)
    if len(missing)==0:
        cover_report.update(exact_upper_bound=380,status='PASS finite 380-clique cover')
        a.out.write_text(json.dumps(cover_report,indent=2)+'\n');print(json.dumps(cover_report),flush=True);return
    if len(N)+16>a.max_pool:
        cover_report.update(status='Incomplete clique cover; unrestricted pair matrix not constructed')
        a.out.write_text(json.dumps(cover_report,indent=2)+'\n');print(json.dumps({k:v for k,v in cover_report.items() if k!='latitude_specs'}),flush=True);return
    E=np.argwhere(np.triu(N.astype(float)@N.astype(float).T>16*np.outer(nd,nd),1))
    rr,cc=np.where(greater3(N@A.T,N@B.T,16*nd[:,None]));E=np.vstack([E,np.c_[rr,n+cc]])
    qr,qc=np.where(np.triu(greater3(A@A.T+3*(B@B.T),A@B.T+B@A.T,16),1));E=np.vstack([E,np.c_[n+qr,n+qc]])
    total=n+len(A)
    M=coo_matrix((np.ones(2*len(E)),(np.repeat(np.arange(len(E)),2),E.ravel())),shape=(len(E),total)).tocsr()
    report=dict(cutoff=a.cutoff,latitude_specs=rows,new_rational_points=len(P),pool_points=total,conflicts=len(E),baseline_sector=380,scope='Finite exact candidate/conflict generation; numerical MILP optimality is not an exact exclusion')
    print(json.dumps({k:v for k,v in report.items() if k!='latitude_specs'}),flush=True)
    np.savez_compressed(a.out.with_suffix('.npz'),rational=N,rational_denominators=nd,quadratic_A=A,quadratic_B=B,edges=E,indices=I)
    res=milp(-np.ones(total),integrality=np.ones(total),bounds=Bounds(0,1),constraints=LinearConstraint(M,-np.inf,1),options={'time_limit':a.seconds,'mip_rel_gap':0,'threads':1})
    report.update(status=int(res.status),message=res.message,seconds=time.monotonic()-start)
    if res.x is not None:
        ids=np.flatnonzero(res.x>.5);assert not np.any(np.isin(E[:,0],ids)&np.isin(E[:,1],ids))
        report.update(selected_count=len(ids),selected_indices=ids.tolist(),mip_gap=float(res.mip_gap),mip_dual_bound=float(res.mip_dual_bound))
    a.out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ['selected_indices','latitude_specs']}),flush=True)


if __name__=='__main__':main()
