"""Exploratory exact-moment models for scaled D_k normal sublattices.

Normal coordinates are h/4.  The dual root lattice forces every h coordinate
to have the same parity.  All cells except the known roots satisfy
max_{i!=j} (|h_i|+|h_j|) <= 8.  Signed permutation averaging is over actual
orientations of the exhibited normal lattice, not ambient shell symmetries.
"""
import itertools as it
import json, math, sys, time
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import numpy as np
from scipy.optimize import linprog, milp, Bounds, LinearConstraint
from scipy.sparse import coo_matrix

N=6218175600
HERE=Path(__file__).resolve().parent

def orbit(rep):
    for p in sorted(set(it.permutations(rep))):
        nz=[i for i,x in enumerate(p) if x]
        for signs in it.product([-1,1],repeat=len(nz)):
            h=list(p)
            for i,s in zip(nz,signs):h[i]*=s
            yield tuple(h)

def partitions(m,k,lo=0):
    if k==0:
        if m==0:yield ()
        return
    for x in range(lo,m+1):
        for rest in partitions(m-x,k-1,x):yield (x,)+rest

def build(k):
    reps=[]
    for parity in [0,1]:
        for h in it.combinations_with_replacement(range(parity,9,2),k):
            if h[-1]+h[-2]<=8:reps.append(h)
    reps=sorted(reps,key=lambda h:(sum(x*x for x in h),h))
    rid={h:i for i,h in enumerate(reps)}
    orbits=[list(orbit(h)) for h in reps]
    parent=list(range(len(reps)))
    def root(x):
        while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
        return x
    def merge(a,b):
        a,b=root(a),root(b)
        if a!=b:parent[max(a,b)]=min(a,b)
    # For a root 2(e_i +/- e_j), dot(v,root)=4 is h_i +/- h_j=8.
    for ir,rr in enumerate(orbits):
        for h in rr:
            for i in range(k):
                for j in range(i):
                    for si,sj in it.product([-1,1],repeat=2):
                        if si*h[i]+sj*h[j]!=8:continue
                        hh=[-x for x in h];hh[i]+=8*si;hh[j]+=8*sj
                        rh=tuple(sorted(abs(x) for x in hh))
                        assert rh in rid
                        merge(ir,rid[rh])
            if k==4:
                hh=(h[0]+h[1]+h[2]+h[3],h[0]+h[1]-h[2]-h[3],h[0]-h[1]+h[2]-h[3],h[0]-h[1]-h[2]+h[3])
                assert all(x%2==0 for x in hh)
                rh=tuple(sorted(abs(x//2) for x in hh));assert rh in rid
                merge(ir,rid[rh])
    group_reps=sorted(set(root(i) for i in range(len(reps))))
    groups=[[h for i,oo in enumerate(orbits) if root(i)==r for h in oo] for r in group_reps]
    gid={h:i for i,gg in enumerate(groups) for h in gg}
    known=[]
    for i in range(k):
        for j in range(i):
            for si,sj in it.product([-1,1],repeat=2):
                h=[0]*k;h[i]=8*si;h[j]=8*sj;known.append(tuple(h))
    rows=[];rhs=[];tags=[]
    for degree in range(6):
        for pp in partitions(degree,k):
            powers=tuple(2*x for x in pp)
            val=lambda h:math.prod(x**p for x,p in zip(h,powers))
            row=[sum(val(h) for h in gg) for gg in groups]
            b=F(N*128**degree*math.prod(math.prod(range(1,2*a,2)) for a in pp),math.prod(72+2*j for j in range(degree)))-sum(val(h) for h in known)
            den=b.denominator;row=[x*den for x in row];bn=b.numerator;g=math.gcd(*row,bn)
            rows.append([x//g for x in row]);rhs.append(bn//g);tags.append(list(powers))
    A=np.asarray(rows,float);b=np.asarray(rhs,float)/N
    colscale=np.array([np.exp(-4.5*sum(x*x for x in gg[0])/16) for gg in groups])
    colscale/=sum(len(gg)*x for gg,x in zip(groups,colscale))
    scales=np.maximum(np.max(np.abs(A*colscale[None,:]),axis=1),np.abs(b))
    As=A*colscale[None,:]/scales[:,None];bs=b/scales
    from scipy.linalg import qr
    _,R,piv=qr(As.T,pivoting=True,mode='economic')
    rank=np.linalg.matrix_rank(As)
    choose=sorted(piv[:rank].tolist())
    print('CONSTRAINT_RANK',k,rank,'of',len(rows),'aug',np.linalg.matrix_rank(np.column_stack([As,bs])),flush=True)
    As=As[choose];bs=bs[choose];scales=scales[choose]
    rows=[rows[i] for i in choose];rhs=[rhs[i] for i in choose];tags=[tags[i] for i in choose]
    base=linprog(np.zeros(len(groups)),A_eq=As,b_eq=bs,bounds=(0,None),method='highs',options={'presolve':False,'threads':1})
    assert base.success,base.message
    print('MODEL',k,'cells',len(gid),'groups',len(groups),'rows',len(rows),'rank',np.linalg.matrix_rank(As),flush=True)
    return dict(k=k,reps=reps,groups=groups,gid=gid,rows=rows,rhs=rhs,tags=tags,As=As,bs=bs,scales=scales,colscale=colscale)

def solve_obj(M,chosen):
    obj=np.zeros(len(M['groups']))
    for h in chosen:obj[M['gid'][tuple(h)]]+=1
    res=linprog(obj*M['colscale'],A_eq=M['As'],b_eq=M['bs'],bounds=(0,None),method='highs',options={'threads':1})
    assert res.success,res.message
    return res.fun*N,res,obj

def serialize_model(M):
    return {key:value for key,value in M.items() if key not in ['As','bs','scales','gid','colscale']}

def search(M,qmax=3,seconds=20):
    k=M['k'];cells=np.array(sorted(M['gid']),dtype=np.int64)
    norm=np.sum(cells*cells,axis=1)
    cells=cells[(norm>0)&(norm<=int(16*qmax))];norm=np.sum(cells*cells,axis=1)
    # Moment-minimal weight for each single orbit cell. Other trial weights can
    # be drawn from feasible worst-case distributions without claiming them.
    weights=[]
    for j in range(len(M['groups'])):
        w=np.zeros(len(M['groups']));w[j]=1
        res=linprog(w*M['colscale'],A_eq=M['As'],b_eq=M['bs'],bounds=(0,None),method='highs',options={'threads':1})
        weights.append(res.fun)
    weights=np.array(weights)
    cellweight=np.array([weights[M['gid'][tuple(h)]] for h in cells])
    reps=[h for h in M['reps'] if 0<sum(x*x for x in h)<=24]
    central=solve_obj(M,[(0,)*k])[0]
    best=dict(value=central,offsets=[(0,)*k],anchor=None,exact_certificate=False)
    print('CENTRAL',k,central,flush=True)
    for anchor in reps:
        a=np.array(anchor);na=a@a;dot=cells@a
        mask=(4*(64-dot)**2<=(128-norm)*(128-na))|(dot>=64)
        pool=cells[mask];pn=norm[mask];w=cellweight[mask]
        m=len(pool)
        dots=pool@pool.T
        good=(4*(64-dots)**2<=(128-pn[:,None])*(128-pn[None,:]))|(dots>=64)
        ii,jj=np.where(np.triu(~good,1));nc=len(ii)
        aa=coo_matrix((np.ones(2*nc),(np.repeat(np.arange(nc),2),np.column_stack([ii,jj]).ravel())),shape=(nc,m)).tocsr()
        exact=False
        for round in range(3):
            if round==0:ww=w
            elif round==1:ww=np.exp(-4.5*pn/16)
            else:ww=np.array([(res.x*M['colscale'])[M['gid'][tuple(h)]] for h in pool])
            if max(ww)<=0:continue
            t=time.time()
            opt=milp(-ww/max(ww),integrality=np.ones(m),bounds=Bounds(np.zeros(m),np.ones(m)),constraints=LinearConstraint(aa,-np.inf,1),options={'time_limit':seconds,'mip_rel_gap':1e-7,'threads':1})
            if opt.x is None:
                print('NO_SOLUTION',k,anchor,'pool',m,'status',opt.status,opt.message,flush=True)
                continue
            selected=pool[opt.x>.5]
            value,res,obj=solve_obj(M,selected)
            print('SEARCH',k,anchor,'pool',m,'conflicts',nc,'round',round,'selected',len(selected),'value',value,'status',opt.status,'seconds',round2(time.time()-t),flush=True)
            if value>best['value']+.01:
                best=dict(value=value,offsets=selected.tolist(),anchor=list(anchor),exact_certificate=False,dual=res.eqlin.marginals.tolist(),objective=obj.tolist())
                (HERE/f'candidate_d{k}.json').write_text(json.dumps(best,indent=2))
                print('BEST',k,value,flush=True)
    return best

def round2(t):return float(f'{t:.2f}')

if __name__=='__main__':
    k=int(sys.argv[1]) if len(sys.argv)>1 else 5
    M=build(k);(HERE/f'model_d{k}.json').write_text(json.dumps(serialize_model(M)))
    print(json.dumps(search(M,qmax=3,seconds=10),indent=2))
