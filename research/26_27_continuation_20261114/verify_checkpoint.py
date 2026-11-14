"""Independent input, reconstruction, norm and conflict-edge checks.

Radical comparisons use integer enclosing intervals, not the search sign code.
This does not certify a new kissing-number lower bound or verify source pairs.
"""
import argparse,hashlib,itertools,json,math
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix

ROOT=Path(__file__).resolve().parent
SOURCES={26:'36732203081894af188e75e2dae69327961b79c7b0781385b1589025b7e0b309',27:'79f8c596a3bdffac1005032af0edaa8790d3c724f146115898237905982e757e'}
BOUNDS={}

def require(x,message):
    if not x:raise AssertionError(message)

def sign_expression(a,b,c=0,d=0,r=3,s=1):
    """Sign of a+b sqrt(r)+c sqrt(s)+d sqrt(rs), by rational intervals."""
    terms={};reduce={1:(1,1),2:(2,1),3:(3,1),4:(1,2),6:(6,1),9:(1,3)}
    for rad,coef in [(1,a),(r,b),(s,c),(r*s,d)]:
        root,factor=reduce[int(rad)];terms[root]=terms.get(root,0)+int(coef)*factor
    terms={k:v for k,v in terms.items() if v}
    if not terms:return 0
    for bits in [16,32,64,128,256]:
        scale=1<<bits;lo=hi=0
        for rad,coef in terms.items():
            key=(rad,bits)
            if key not in BOUNDS:
                low=math.isqrt(rad*scale*scale);BOUNDS[key]=(low,low if low*low==rad*scale*scale else low+1)
            lower,upper=BOUNDS[key];lo+=coef*(lower if coef>0 else upper);hi+=coef*(upper if coef>0 else lower)
        if lo>0:return 1
        if hi<0:return -1
    raise AssertionError('Undecided exact radical comparison')

def source(dim):
    n=201180 if dim==26 else 209862;path=ROOT/f'data/PackingStar/25D-31D/new_bounds_20261111/D{dim}_{n}_coordinates.npz'
    require(hashlib.sha256(path.read_bytes()).hexdigest()==SOURCES[dim],'Source hash')
    z=np.load(path,allow_pickle=False)
    if dim==26:return 2*np.r_[z['rational_points'],z['quadratic_A']],2*np.r_[np.zeros_like(z['rational_points']),z['quadratic_B']]
    return 2*z['points'],np.zeros_like(z['points'])

def norms(A,B,rad,den):
    require(A.shape==B.shape and A.shape[1]<=27,'Coordinate shape');require(set(map(int,rad))<={1,2,3},'Fields');require(np.all((den>0)&(den<=32767)),'Denominators');require(np.max(np.abs(A.astype(np.int64)))<=32767 and np.max(np.abs(B.astype(np.int64)))<=32767,'Coefficient bound')
    aa=A.astype(np.int64);bb=B.astype(np.int64)
    require(np.all(np.sum(aa*bb,axis=1)==0),'Irrational norm coefficient');require(np.all(np.sum(aa*aa+rad[:,None]*bb*bb,axis=1)==128*den.astype(np.int64)**2),'Rational norm coefficient')

def baseline(dim):
    z=np.load(ROOT/f'data/cycle{dim}/data.npz',allow_pickle=False);L=z['leech'];T=z['tail_palette'];flat=z['state'][z['initial']];ci=flat//len(T);ph=flat%len(T);H=z['H'].reshape(-1,4,24)[ci];C=z['offset'].reshape(-1,4)[ci]^ph[:,None];P=z['parents'][ci];look={v.tobytes():i for i,v in enumerate(L)};anti=np.array([look[(-v).tobytes()] for v in L]);keep=~z['blocked'].copy();keep[P.ravel()]=False;keep[anti[P.ravel()]]=False
    LC=np.c_[2*L[keep],np.zeros((int(keep.sum()),dim-24),np.int16)];middle=np.array([np.r_[s*h,t*T[c]] for hs,cs in zip(H,C) for h,c in zip(hs,cs) for s,t in itertools.product([-1,1],repeat=2)],np.int16);A=np.r_[LC,z['fixed_A'],middle];B=np.r_[np.zeros_like(LC),z['fixed_B'],np.zeros_like(middle)];SA,SB=source(dim);keys={np.r_[a,b].tobytes() for a,b in zip(A,B)}
    require(len(keys)==len(A)==len(SA),'Baseline count/distinctness');require(keys=={np.r_[a,b].tobytes() for a,b in zip(SA,SB)},'Baseline exact set reconstruction');return len(A)

def one_missing(dim,wrong_cost=False,wrong_parent=False):
    z=np.load(ROOT/f'data/one_missing{dim}.npz',allow_pickle=False);A=z['candidate_A'];P=z['parent_indices'].copy();reported=z['orbit_blocker_counts'].copy();L=np.load(ROOT/'data/leech.npy',allow_pickle=False);g=np.load(ROOT/f'data/parentgraph{dim}/graph.npz',allow_pickle=False);O,OB=source(dim)
    if wrong_cost:reported[0]+=1
    if wrong_parent:P[0,0]=P[0,1]
    norms(A,np.zeros_like(A),np.ones(len(A),int),np.ones(len(A),int));H=A[:,:24];T=A[:,24:]
    for h,ps in zip(H,P):
        dots=L[ps].astype(np.int64)@h.astype(np.int64);require(np.all(np.abs(dots)==48),'Leech parent products');require(np.array_equal(np.sum(L[ps]*np.sign(dots)[:,None],axis=0),h),'Parent sum');require(np.sum(g['source_for_leech'][ps]<0)==1,'One absent parent')
    oldH=g['bodies'];positive=oldH[np.arange(len(oldH)),np.argmax(oldH!=0,axis=1)]>0;oldH=oldH[positive];oldT=g['tails'][positive];mid=np.sum(np.abs(H.astype(np.int64)@oldH.T)+np.abs(T.astype(np.int64)@oldT.T)>64,axis=1);core=np.sum(O[:,:24].astype(np.int64)**2+3*OB[:,:24].astype(np.int64)**2,axis=1);groups={};quads={}
    for i in np.flatnonzero((core==64)|np.any(OB,axis=1)):
        q=OB[i,:24].any();h=(OB[i,:24] if q else O[i,:24]).copy();t=O[i,24:].copy()
        if h[np.flatnonzero(h)[0]]<0:h=-h
        if t[np.flatnonzero(t)[0]]<0:t=-t
        (quads if q else groups)[np.r_[h,t].tobytes()]=(h,t)
    HA=np.array([v[0] for v in groups.values()],np.int64);HT=np.array([v[1] for v in groups.values()],np.int64);count=mid+np.sum(np.abs(H.astype(np.int64)@HA.T)+np.abs(T.astype(np.int64)@HT.T)>64,axis=1)
    for qb,qt in quads.values():
        p=np.abs(T.astype(np.int64)@qt)-64;q=np.abs(H.astype(np.int64)@qb)
        count+=np.array([sign_expression(a,b)>0 for a,b in zip(p,q)])
    require(np.array_equal(count,reported),'Orbit blocker count');return len(A)

def graph_case():
    z=np.load(ROOT/'data/point_case27/nodes.npz',allow_pickle=False);A=z['A'];B=z['B'];rad=z['radicands'];den=z['denominators'];norms(A,B,rad,den)
    with open(ROOT/'data/point_case27/graph.bin','rb') as f:
        n,m=map(int,np.fromfile(f,'<u4',2));ne=int(np.fromfile(f,'<u8',1)[0]);initial=np.fromfile(f,'<u4',m);ptr=np.fromfile(f,'<u8',n+1);idx=np.fromfile(f,'<u4',ne);require(not f.read(1),'Graph trailing bytes')
    require(n==len(A) and np.array_equal(initial,np.arange(m)),'Graph header');G=csr_matrix((np.ones(ne,np.int8),idx.astype(np.int32),ptr.astype(np.int64)),shape=(n,n));require((G!=G.T).nnz==0 and not G.diagonal().any(),'Graph symmetry');require(G[:m,:m].nnz==0,'Incumbent independence');SA,SB=source(27);known={(1,1,a.tobytes(),b.tobytes()) for a,b in zip(SA,SB)};keys=[(int(r),int(d),a.tobytes(),b.tobytes()) for a,b,r,d in zip(A,B,rad,den)];require(len(set(keys))==n,'Node distinctness');require(all(k in known for k in keys[:m]),'Incumbent source membership');require(all(k not in known for k in keys[m:]),'Candidate novelty');ii,jj=G.nonzero();up=ii<jj;ii=ii[up];jj=jj[up];aa=A.astype(np.int64);bb=B.astype(np.int64)
    for at in range(0,len(ii),4096):
        u=ii[at:at+4096];v=jj[at:at+4096];p=np.einsum('ij,ij->i',aa[u],aa[v])-64*den[u].astype(np.int64)*den[v];q=np.einsum('ij,ij->i',bb[u],aa[v]);c=np.einsum('ij,ij->i',aa[u],bb[v]);e=np.einsum('ij,ij->i',bb[u],bb[v])
        for a,b,c,d,r,s in zip(p,q,c,e,rad[u],rad[v]):require(sign_expression(a,b,c,d,int(r),int(s))>0,'False conflict edge')
    return dict(nodes=n,old=m,candidates=n-m,verified_undirected_edges=len(ii)),(A,B,rad,den)

def module604(duplicate=False):
    raw=(ROOT/'data/module604_11.json').read_bytes();require(hashlib.sha256(raw).hexdigest()=='0bde9ca2c434d7beb5af63a64a291498e8c264d9a5718b738e92dd70b9ba7761','604 input hash');V=np.array(json.loads(raw)['vectors'],np.int64);P=V[:,::2];Q=V[:,1::2]
    if duplicate:P[1]=P[0];Q[1]=Q[0]
    require(np.all(np.sum(P*P+2*Q*Q,axis=1)==36) and not np.any(np.sum(P*Q,axis=1)),'604 norms');G=P@P.T+2*Q@Q.T;H=P@Q.T+Q@P.T
    for i in range(len(P)):
        for j in range(i):require(sign_expression(int(G[i,j])-18,int(H[i,j]),r=2)<=0,'604 pair')
    return len(P)*(len(P)-1)//2

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--controls',action='store_true');ap.add_argument('--out',type=Path,default=ROOT/'verification.json');a=ap.parse_args();manifest=json.loads((ROOT/'manifest.json').read_text())
    for rel,digest in manifest['files'].items():require(hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==digest,f'Artifact hash: {rel}')
    L=np.load(ROOT/'data/leech.npy',allow_pickle=False);require(L.shape==(196560,24) and np.all(np.sum(L.astype(np.int64)**2,axis=1)==32),'Leech seed norms');counts={d:baseline(d) for d in [26,27]};one={d:one_missing(d) for d in [26,27]};graph,arrays=graph_case();pairs=module604();controls=0
    if a.controls:
        def reject(fn):
            nonlocal controls
            try:fn()
            except AssertionError:controls+=1;return
            raise AssertionError('Corruption was accepted')
        AA,BB,rr,dd=arrays;bad=AA.copy();bad[0,0]+=1;reject(lambda:norms(bad,BB,rr,dd));reject(lambda:one_missing(26,wrong_cost=True));reject(lambda:one_missing(26,wrong_parent=True));reject(lambda:module604(duplicate=True))
        p=int(AA[0].astype(np.int64)@AA[1]);reject(lambda:require(sign_expression(p-64*int(dd[0])*int(dd[1]),0,r=1)>0,'Compatible incumbent pair'))
    report=dict(new_lower_bound=False,source_all_pairs_checked=False,hashed_files=len(manifest['files']),baseline_reconstruction=counts,one_missing_states=one,point_case27=graph,module604_verified_pairs=pairs,rejected_corruptions=controls);a.out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
