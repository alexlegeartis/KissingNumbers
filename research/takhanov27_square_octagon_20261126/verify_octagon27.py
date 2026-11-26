import numpy as np, itertools, sys, collections, hashlib
from pathlib import Path
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
D=Path(sys.argv[1]).resolve()
if (D/'common.py').is_file() and (D/'data/d27_certificate.npz').is_file():
    # Native source layout: k-nic/KissingNumbers@77f0860b/version1/certificates.
    cf=(D/'common.py').read_bytes()
    assert hashlib.sha1(b'blob '+str(len(cf)).encode()+bytes([0])+cf).hexdigest()=='d1ae49a0f3ab6e9fb9906a98592ab9dbc6d9d7d5'
    sys.path.insert(0,str(D))
    from common import leech_shell
    fp=D/'data/d27_certificate.npz'
    assert hashlib.sha256(fp.read_bytes()).hexdigest()=='8b704bb49503e3a4f6ad6b567690a7a2a9abfb30600f453313eb0510e6d29a94'
    f=np.load(fp,allow_pickle=False);L=leech_shell().astype('<i2')
    keep=np.ones(len(L),dtype=bool);keep[f['deleted_leech_indices']]=False
    Q=np.vstack([np.pad(L[keep],((0,0),(0,3))),f['added_integer'].astype('<i2')])
else:
    for filename,digest in {'d27_certificate.npz':'8b704bb49503e3a4f6ad6b567690a7a2a9abfb30600f453313eb0510e6d29a94','L.npy':'68b3f12fa1d7024c6e3fafd91b3ea616526ccb1299b76c3b6e6e7d53d89ff302','Q.npy':'9674e0dc67eebf83be574c7da07292cd98a8a6c533b4d786591bb08390165d3b'}.items():
        assert hashlib.sha256((D/filename).read_bytes()).hexdigest()==digest
    f=np.load(D/'d27_certificate.npz');L=np.load(D/'L.npy');Q=np.load(D/'Q.npy')
assert np.array_equal(2*f['added_integer'][f['midpoint_rows'],:24],L[f['midpoint_endpoints']].sum(1))
assert len(L)==196560 and L.shape[1]==24
assert len(Q)==210108 and np.array_equal(Q[-26796:],f['added_integer']) and np.all(np.sum(Q.astype('i4')**2,1)==128)
ids=collections.defaultdict(list)
for i,(a,b) in zip(f['midpoint_rows'],f['midpoint_endpoints']):
    ids[tuple(sorted((int(a),int(b))))].append(i)
edges=list(ids);inc=collections.defaultdict(list)
for e,(a,b) in enumerate(edges):inc[a].append(e);inc[b].append(e)
cs=[];seen=set()
for r in sorted(inc):
    if r in seen:continue
    u=r;last=None;c=[]
    while u not in seen:
        seen.add(u);c.append(u);i=next(i for i in inc[u] if i!=last);u=sum(edges[i])-u;last=i
    if len(c)==8:cs.append(c)
X=Q[-26796:];E0=len(Q)-len(X); exterior=set(map(int,f['deleted_leech_indices']));global_id={k:i for i,k in enumerate(j for j in range(len(L)) if j not in exterior)}
def key(a):return a.astype('<i2').tobytes()
lookup={key(x):i for i,x in enumerate(L)}
T=np.array([(a,b,c) for a,b,c in itertools.product([-4,0,4],repeat=3) if a*a+b*b+c*c==32 and next(t for t in (a,b,c) if t)>0],dtype='i2')
MATRIX=[a for a in itertools.product(list(itertools.combinations(range(4),2)),repeat=4) if all(sum(j in r for r in a)==2 for j in range(4))]
assert len(MATRIX)==90
whole=[];bag={};local=[];histo=collections.Counter()
for ci,c in enumerate(cs):
    R=L[c]
    if (R[4]==-R[0]).all():continue
    grams=np.array([[128 if i==j else 64 if (i-j)%8 in [1,7] else 0 for j in range(8)] for i in range(8)])
    assert np.array_equal(R@R.T,grams)
    U=np.zeros((4,4,24),dtype='i2');uix=np.zeros((4,4),int)
    for i,j in itertools.product(range(4),repeat=2):
        start=2*i;end=2*j+1
        U[i,j]=sum((-1)**k*R[(k+start)%8] for k in range((end-start)%8))
        uix[i,j]=lookup[key(U[i,j])]
    drop=sum(int(x) in exterior for x in uix.flatten())
    bodies=[];ends=[];bix={}
    for hv,l in itertools.product(range(2),range(4)):
        for j,k in itertools.combinations(range(4),2):
            a,b=((l,j),(l,k)) if hv==0 else ((j,l),(k,l))
            assert U[a]@U[b]==64 and np.all((U[a]+U[b])%2==0)
            bix[hv,l,j,k]=len(bodies);bodies.append((U[a]+U[b])//2);ends.append([int(uix[a]),int(uix[b])])
    B=np.array(bodies,dtype='i2');assert np.all(np.sum(B.astype('i4')**2,1)==96) and np.all(np.sum(T.astype('i4')**2,1)==32)
    P=(B.astype('f')@X[:,:24].astype('f').T).astype('i2')
    Ft=(T.astype('f')@X[:,24:].astype('f').T).astype('i2')
    C={};signC={}
    def reg(point,guards):
        p=key(point);v=frozenset(guards)
        if p in bag: assert bag[p]==v
        else:bag[p]=v
    for b,t in itertools.product(range(48),range(6)):
        g0=set(global_id[s] for s in ends[b] if s not in exterior)
        gs=[]
        for s in [-1,1]:
            guards=g0|set(map(int,np.flatnonzero(P[b]+s*Ft[t]>64)+E0))
            reg(np.r_[B[b],s*T[t]],guards);gs.append(guards)
        C[b,t]=gs[0]|gs[1];signC[b,t]=gs
    best=(-999999,None);bestpartial=0;hdef=collections.Counter()
    Prod=(U.reshape(16,24).astype('f')@X[:,:24].astype('f').T).astype('i2')
    refunds=[set(map(int,np.flatnonzero(Prod[j]>64)+E0)) for j,ix in enumerate(uix.flat) if int(ix) in exterior]
    for axis in range(3):
        t0,t1=np.where(T[:,axis]==0)[0]
        for mat in MATRIX:
            block=[]
            for hv,l in itertools.product(range(2),range(4)):
                on=set(mat[l]) if hv==0 else set(i for i in range(4) if l in mat[i])
                off=set(range(4))-on
                for state,jk in ((1,on),(0,off)):
                    kk=bix[hv,l,*sorted(jk)]
                    t=(t0,t1)[state^hv];block.append((kk,t))
            rows=np.array([np.r_[B[b],s*T[t]] for b,t in block for s in [-1,1]])
            gram=rows@rows.T
            assert np.all(np.diag(gram)==128) and np.all(np.triu(gram,1)<=64)
            kill=set.union(*(C[p] for p in block));gain=32-len(kill)
            own=[z for p in block for z in signC[p]]+refunds
            match={}
            def aug(i,seen):
                for q in own[i]:
                    if q in seen:continue
                    seen.add(q)
                    if q not in match or aug(match[q],seen):match[q]=i;return True
                return False
            for ii in sorted(range(len(own)),key=lambda i:len(own[i])):aug(ii,set())
            hdef[len(own)-len(match)]+=1
            if len(match)<len(own):
                points=np.array([np.r_[B[b],s*T[t]] for b,t in block for s in [-1,1]]+[np.r_[r,np.zeros(3,dtype='i2')] for r,ix in zip(U.reshape(16,24),uix.flat) if int(ix) in exterior])
                G=points@points.T>64;np.fill_diagonal(G,False)
                memo=set()
                def solve(pool):
                    if pool in memo:return None
                    memo.add(pool);mm={}
                    def reach(i,seen):
                        for q in own[i]:
                            if q in seen:continue
                            seen.add(q)
                            if q not in mm or reach(mm[q],seen):mm[q]=i;return True
                        return False
                    for i in sorted(pool,key=lambda x:len(own[x])):reach(i,set())
                    left=set(pool)-set(mm.values())
                    if not left:return None
                    todo=list(left);right=set()
                    while todo:
                        j=todo.pop()
                        for q in own[j]-right:
                            right.add(q)
                            if q in mm and mm[q] not in left:left.add(mm[q]);todo.append(mm[q])
                    assert len(left)>len(right)
                    l=sorted(left);gg=G[np.ix_(l,l)]
                    if not gg.any():
                        return (l,sorted(right))
                    z=int(np.argmax(gg.sum(1)));az=l[z];bz=l[next(i for i in np.argsort(-gg.sum(1)) if gg[z,i])]
                    return solve(pool-{az}) or solve(pool-{bz})
                witness=solve(frozenset(range(len(own))))
                if witness is not None:
                    ind,deleted=witness;fresh=points[ind]
                    bb=fresh.astype('f')@Q.astype('f').T
                    assert set(map(int,np.flatnonzero(np.any(bb>64,axis=0))))==set(deleted)
                    assert np.all((np.triu(fresh@fresh.T,1))<=64)
                    np.savez('/tmp/lower48_oct_win.npz',added=fresh,delete=deleted)
                    raise RuntimeError('WIN '+str((ci,len(ind),len(deleted))))
            if gain>best[0]:best=(gain,[axis,mat,len(kill),len(set(k for k in kill if k<E0))])
    # Include original K4,4 cells too, so the global matching covers arbitrary refunds.
    Prod=(U.reshape(16,24).astype('f')@X[:,:24].astype('f').T).astype('i2')
    for j,ix in enumerate(uix.flat):
        reg(np.r_[U.reshape(16,24)[j],np.zeros(3,'i2')],set(map(int,np.flatnonzero(Prod[j]>64)+E0))|({global_id[int(ix)]} if int(ix) in global_id else set()))
    assert dict(hdef)=={0:249,1:8,2:12,8:1};print('chart',ci,'all 270 incl refunds excluded',flush=True)
    whole.append(best[0]);histo[drop]+=1;local.append([ci,best])
def Hall(bag):
    common={k:i for i,k in enumerate(sorted(set().union(*bag.values())))};pos=[0];ids=[]
    for v in bag.values():ids.extend(common[k] for k in sorted(v));pos.append(len(ids))
    G=csr_matrix((np.ones(len(ids),dtype='b'),ids,pos),shape=(len(bag),len(common)));match=maximum_bipartite_matching(G,perm_type='column')
    return [len(bag),len(common),len(ids),int(sum(match>=0))]
print('checked',len(whole),'A7-carriers; deleted-cell histogram',dict(histo),'best per-carrier 16-body full matrix=',collections.Counter(whole),'best detail',max(local,key=lambda x:x[1][0]))
print('all donors and original cells across all carriers [N,guards,edges,matched]',Hall(bag))
