#!/usr/bin/env python3
"""Independent reconstruction of Li's 12460-point 19D configuration.

Source: Section 4 and Table 5, https://hexagonmath.org/pdf/2610.00248v1.
Every geometric compatibility decision is exact in Z[sqrt(6)].
"""
import itertools,json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent;DATA=HERE/'data'
GENERATORS=[(1,8,9,12,16,17,18,19),(2,10,11,14,15,17,18,20),(3,7,9,13,15,16,17,19),
 (4,7,8,10,12,15,16,19),(5,10,12,13,15,16,17,18),(6,7,8,9,10,13,16,18),
 (1,4,7,9),(1,5,6,18,20),(1,3,12,15,20),(1,10,13,19,20),
 (1,3,5,6,7,13,14,15,18,20),(2,4,6,7,8,13,14,16,17,18)]
REPS_AXIS=[1134,2453,4337,5590,6500,7434,8363,9413,10359,11600,12378,15592,17037,18346,19224,20342,
 21301,24455,26184,27571,29399,30393,31243,32556,66280,69210,70169,71287,72389,73698,74307,75620,
 76758,77752,79580,80679,82021,82955,84153,85406,87290,88321,90503,93493,94271,95512,96682,97732]
REPS_DIAG=[986,1972,2822,3617,4907,5644,6846,7888,9073,10015,10980,12227,13184,13991,15196,16178,
 16414,17520,18626,19941,20902,21633,22579,23645,25084,26002,26729,27982,28740,30051,30872,31990,
 65842,66581,67751,68809,70083,71085,71967,72760,73936,75255,75788,76898,77857,78927,80308,81043,
 82678,83921,84835,85773,86862,87840,88978,89781,90717,92026,92801,93935,95205,96139,96880,98135]


def span(basis):
    S=[0]
    for b in basis:S += [w^b for w in S]
    return S


def nullspace(rows,n):
    rows=list(rows);pivots=[];r=0
    for p in range(n):
        k=next((i for i in range(r,len(rows)) if rows[i]>>p&1),None)
        if k is None:continue
        rows[r],rows[k]=rows[k],rows[r]
        for i in range(len(rows)):
            if i!=r and rows[i]>>p&1:rows[i]^=rows[r]
        pivots.append(p);r+=1
    out=[]
    for p in range(n):
        if p not in pivots:
            out.append((1<<p)|sum(((rows[i]>>p)&1)<<q for i,q in enumerate(pivots)))
    return out


def le0(a,b):
    a,b=np.broadcast_arrays(np.asarray(a,dtype=np.int64),np.asarray(b,dtype=np.int64))
    out=(b==0)&(a<=0)
    out|=(b<0)&(a<=0)
    for mask,ge in (((b>0)&(a<=0),True),((b<0)&(a>0),False)):
        aa=a[mask];bb=b[mask]
        if len(aa):
            if max(map(abs,map(int,aa)))**2>=2**63 or 6*max(map(abs,map(int,bb)))**2>=2**63:
                aa=aa.astype(object);bb=bb.astype(object)
            out[mask]=(aa*aa>=6*bb*bb) if ge else (aa*aa<=6*bb*bb)
    return out


def verify_blocks(blocks):
    # Each block is (r,s,Q) representing (r+s sqrt(6))/sqrt(Q).
    count=0;bad=0;contacts=0
    for i,(R,S,Q) in enumerate(blocks):
        for j,(U,V,P) in enumerate(blocks):
            if j<i:continue
            for a in range(0,len(R),256):
                rr=R[a:a+256].astype(float);ss=S[a:a+256].astype(float)
                A=np.rint(rr@U.astype(float).T+6*ss@V.astype(float).T).astype(np.int64)
                B=np.rint(rr@V.astype(float).T+ss@U.astype(float).T).astype(np.int64)
                product=int(Q)*int(P)
                C=4*A*A+24*B*B-product;D=8*A*B
                valid=le0(A,B)|le0(C,D)
                if i==j:
                    for t in range(len(A)):valid[t,a+t]=True
                bad+=int(np.sum(~valid))
                contacts+=int(np.sum((C==0)&(D==0)))
                count+=valid.size
    assert bad==0,bad
    return dict(checked_matrix_entries=count,violations=bad,rational_contacts_counted_with_symmetry=contacts)


def prepare():
    DATA.mkdir(exist_ok=True)
    gens=[sum(1<<(p-1) for p in row) for row in GENERATORS]
    E=span(gens);assert len(set(E))==4096
    dual=span(nullspace(gens,20));assert len(dual)==256
    octads=[w for w in dual if w.bit_count()==8];assert len(octads)==130
    T=(0,3,6,8);Tm=sum(1<<p for p in T);H=[p for p in range(20) if p not in T];Hm=((1<<20)-1)^Tm
    Rmap=np.array([[1,1,-1,-1],[1,-1,1,-1],[1,-1,-1,1]],dtype=np.int64)
    tetrads=sorted(w for w in E if w.bit_count()==4);assert len(tetrads)==5 and Tm in tetrads
    U0=[w for w in span(tetrads) if sum(bool(w&t) for t in tetrads)%2==0]
    assert len(U0)==16
    W=sorted({r^u for r in REPS_AXIS+REPS_DIAG for u in U0});assert len(W)==1792
    Ws=set(W);assert Ws<=set(E) and all(w^Tm not in Ws for w in W)
    rows=[]
    for i,j in itertools.combinations(range(20),2):
        for si,sj in itertools.product((-2,2),repeat=2):
            v=[0]*20;v[i]=si;v[j]=sj;rows.append(v)
    for word in octads:
        positions=[i for i in range(20) if word>>i&1]
        for bits in range(256):
            if bits.bit_count()%2:
                v=[0]*20
                for j,p in enumerate(positions):v[p]=-1 if bits>>j&1 else 1
                rows.append(v)
    rows=np.asarray(rows,dtype=np.int64);section=rows[rows[:,T].sum(1)==0]
    tail=section[:,T]@Rmap.T;assert np.all(tail%2==0)
    base=np.concatenate([section[:,H],tail//2],axis=1).astype(np.int16)
    assert base.shape==(10668,19) and np.all((base.astype(np.int64)**2).sum(1)==8)
    low=[w for w in E if (w&Hm).bit_count()==4 and (w&Tm).bit_count()==1];assert len(low)==64
    axis=[];diagR=[];diagS=[];marks=[];axiswords=[];diagwords=[]
    for word in W:
        q=np.asarray([1-2*((word>>p)&1) for p in range(20)],dtype=np.int64)
        tail=Rmap@q[list(T)];assert np.all(tail%2==0);tail//=2
        kind=(word&Tm).bit_count()
        if kind==2:
            differences=[d for d in low if word^d in Ws]
            occupied=0
            for d in differences:occupied|=d
            free=[p for p in H if not(occupied>>p&1)]
            assert len(differences)==24 and len(free)==1
            mark=H.index(free[0]);v=np.concatenate([500*q[H],500*tail]);v[mark]=373*q[H[mark]]
            axis.append(v);marks.append(mark);axiswords.append(word)
        else:
            assert kind in (1,3) and np.all(np.abs(tail)==1)
            r=np.concatenate([3*q[H],np.zeros(3,dtype=int)]);s=np.concatenate([np.zeros(16,dtype=int),tail])
            diagR.append(r);diagS.append(s);diagwords.append(word)
    axis=np.asarray(axis,dtype=np.int16);diagR=np.asarray(diagR,dtype=np.int16);diagS=np.asarray(diagS,dtype=np.int16)
    assert len(axis)==768 and len(diagR)==1024
    assert np.all((axis.astype(np.int64)**2).sum(1)==4889129)
    assert np.all((diagR.astype(np.int64)**2+6*diagS.astype(np.int64)**2).sum(1)==162)
    np.savez_compressed(DATA/'baseline19.npz',base=base,axis=axis,diagR=diagR,diagS=diagS,
                        E=E,octads=octads,W=W,axis_words=axiswords,diag_words=diagwords,marks=marks,
                        head_positions=H,tetrad_positions=T,tetrads=tetrads)
    blocks=[(base,np.zeros_like(base),8),(axis,np.zeros_like(axis),4889129),(diagR,diagS,162)]
    report=dict(dimension=19,baseline=12460,base=10668,axis=768,diagonal=1024,
                exact_field='Z[sqrt(6)]',**verify_blocks(blocks),source='https://hexagonmath.org/pdf/2610.00248v1')
    (DATA/'preparation19.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':prepare()
