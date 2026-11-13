"""Exact 438-point E8/D3 module and current-source 27D union tests.

240 E8 roots + 3*4*16 caps + 6 pure tails =438. For each D3 tail plane,
use a different16-point head cross-polytope; distinct frames have all absolute
cross products<=8. Coordinates are in Q(sqrt2), with exact norm32.
"""
from pathlib import Path
from itertools import combinations,product
import argparse,json,time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import maximum_flow,breadth_first_order


def greater2(p,q,lim=16):
    p=np.asarray(p,np.int64);q=np.asarray(q,np.int64);r=lim-p
    return np.where(q>=0,(r<0)|(2*q*q>r*r),(r<0)&(r*r>2*q*q))


def module():
    H4=np.asarray([[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]],np.int16)
    F=[np.vstack([4*np.eye(8,dtype=np.int16),-4*np.eye(8,dtype=np.int16)])]
    for blocks in [[[0,1,2,3],[4,5,6,7]],[[0,1,4,5],[2,3,6,7]]]:
        r=[]
        for block in blocks:
            for row in H4:
                h=np.zeros(8,np.int16);h[block]=2*row;r.extend([h,-h])
        F.append(np.asarray(r))
    for i in range(3):
        g=F[i].astype(int)@F[i].astype(int).T;assert np.all(np.diag(g)==16);np.fill_diagonal(g,0);assert g.max()<=0
        for j in range(i):assert np.abs(F[i].astype(int)@F[j].astype(int).T).max()<=8
    A=[];B=[]
    for i,j in combinations(range(8),2):
        for s,t in product([-4,4],repeat=2):
            r=np.zeros(11,np.int16);r[i]=s;r[j]=t;A.append(r);B.append(np.zeros(11,np.int16))
    for s in product([-2,2],repeat=8):
        if sum(x<0 for x in s)%2:continue
        r=np.zeros(11,np.int16);r[:8]=s;A.append(r);B.append(np.zeros(11,np.int16))
    for frame,(i,j) in zip(F,combinations(range(3),2)):
        for s,t in product([-2,2],repeat=2):
            for h in frame:
                r=np.zeros(11,np.int16);r[:8]=h;b=np.zeros(11,np.int16);b[8+i]=s;b[8+j]=t;A.append(r);B.append(b)
    for i in range(3):
        for s in [-4,4]:
            b=np.zeros(11,np.int16);b[8+i]=s;A.append(np.zeros(11,np.int16));B.append(b)
    A=np.asarray(A);B=np.asarray(B);assert len(A)==438
    assert np.all(np.sum(A.astype(int)**2+2*B.astype(int)**2,axis=1)==32) and not np.any(np.sum(A.astype(int)*B.astype(int),axis=1))
    p=A.astype(int)@A.astype(int).T+2*(B.astype(int)@B.astype(int).T);q=A.astype(int)@B.astype(int).T+B.astype(int)@A.astype(int).T;np.fill_diagonal(p,-32)
    assert not greater2(p,q).any()
    return A,B


def orient(A,B,trial,rng,rotate=False):
    A=A.copy();B=B.copy()
    if trial:
        perm=rng.permutation(8);s=rng.choice([-1,1],8);s[-1]=np.prod(s[:-1]);A[:,:8]=A[:,perm]*s
        for _ in range(1+trial%5):
            s=rng.choice([-1,1],8);s[-1]=np.prod(s[:-1]);dots=A[:,:8].astype(int)@s;assert np.all(dots%4==0)
            A[:,:8]-=(dots//4)[:,None]*s
        perm=rng.permutation(3);s=rng.choice([-1,1],3);B[:,8:]=B[:,8+perm]*s
    den=1
    if rotate:
        q=np.asarray([0,1,1,1]) if trial==0 else rng.integers(-3,4,size=4)
        if not np.any(q):q[0]=1
        w,x,y,z=map(int,q);den=w*w+x*x+y*y+z*z
        R=np.asarray([[w*w+x*x-y*y-z*z,2*(x*y-w*z),2*(x*z+w*y)],
                      [2*(x*y+w*z),w*w-x*x+y*y-z*z,2*(y*z-w*x)],
                      [2*(x*z-w*y),2*(y*z+w*x),w*w-x*x-y*y+z*z]],np.int64)
        assert np.array_equal(R.T@R,den*den*np.eye(3,dtype=np.int64))
        A=A.astype(np.int64)*den;B=B.astype(np.int64);B[:,8:]=B[:,8:]@R
    assert np.all(np.sum(A.astype(np.int64)**2+2*B.astype(np.int64)**2,axis=1)==32*den*den)
    return A,B,den


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--trials',type=int,default=256)
    ap.add_argument('--out',type=Path,required=True);ap.add_argument('--rotate',action='store_true');a=ap.parse_args();start=time.monotonic();root=Path(__file__).resolve().parents[3]
    src=root.parent/'PackingStar/25D-31D/new_bounds_20261111/D27_209862_coordinates.npz';z=np.load(src);X=z['points'].astype(np.int16);assert int(z['coordinate_denominator'])==1
    J=np.r_[[0,2,4,5,6,10,11,23],np.arange(24,27)];outside=np.any(X[:,[i for i in range(27) if i not in J]],axis=1)
    profiles,inverse,weights=np.unique(X[:,J],axis=0,return_inverse=True,return_counts=True)
    ext_weights=np.bincount(inverse[outside],minlength=len(profiles));A,B=module();rng=np.random.default_rng(20261113);rows=[]
    np.savez_compressed(a.out.with_name('lambda11_438_module.npz'),rational=A,radical=B,radicand=2,squared_norm=32,inner_product_limit=16)
    for trial in range(a.trials):
        C,D,den=orient(A,B,trial,rng,a.rotate);bad=greater2(C.astype(np.int32)@profiles.astype(np.int32).T,D.astype(np.int32)@profiles.astype(np.int32).T,16*den)
        rr,cc=np.where(bad);m=len(C);n=len(profiles);s=m+n;t=s+1
        I=np.r_[np.full(m,s),rr,m+np.arange(n)];K=np.r_[np.arange(m),m+cc,np.full(n,t)];V=np.r_[np.ones(m,np.int32),np.full(len(rr),m+1,np.int32),weights.astype(np.int32)]
        G=coo_matrix((V,(I,K)),shape=(t+1,t+1)).tocsr();flow=maximum_flow(G,s,t);gain=m-int(flow.flow_value);ext=int(ext_weights[np.any(bad,axis=0)].sum())
        row=dict(trial=trial,module438=True,denominator=den,old_module=int((~outside).sum()),external_blockers=ext,whole_replacement_gain=16-ext,exact_union_gain=gain,seconds=time.monotonic()-start)
        if gain>0:
            residual=G-flow.flow;residual.data[residual.data<=0]=0;residual.eliminate_zeros();reach=breadth_first_order(residual,s,directed=True,return_predecessors=False)
            sel=reach[reach<m];deleted=reach[(reach>=m)&(reach<m+n)]-m;keep=~np.isin(inverse,deleted)
            AA=np.zeros((len(sel),27),np.int32);BB=np.zeros_like(AA);AA[:,J]=C[sel];BB[:,J]=D[sel]
            np.savez_compressed(a.out.with_name(a.out.stem+f'_gain{gain}_trial{trial}.npz'),rational=np.vstack([den*X[keep].astype(np.int32),AA]),radical=np.vstack([np.zeros_like(X[keep]),BB]),radicand=2,squared_norm=32*den*den,inner_product_limit=16*den*den,coordinate_denominator=den)
            row.update(removed=int((~keep).sum()),added=len(sel))
        rows.append(row)
        if trial%20==0 or gain:print(json.dumps(row),flush=True)
        a.out.write_text(json.dumps(dict(source=str(src),module_points=438,baseline_points=len(X),best_gain=max(r['exact_union_gain'] for r in rows),best_external_cost=min(r['external_blockers'] for r in rows),scope='438 module exact; union mincuts exact; complete final source replay remains required for a positive result',trials=rows),indent=2)+'\n')


if __name__=='__main__':main()
