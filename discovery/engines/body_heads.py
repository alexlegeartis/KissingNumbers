"""Finite body/head product search with exact gates on each emitted exchange.

Body coordinates and complete-family upper envelopes are supplied inputs.
The result is a proposal, never an independently certified full configuration.
"""
from collections import defaultdict
from itertools import combinations
from math import isfinite
from time import monotonic

import numpy as np

from ..core import Artifact, Verification, check_params, integer, digest


def _matrix(raw, width, name, context, empty=False):
    if not isinstance(raw,list) or (not raw and not empty) or any(
            not isinstance(row,list) or len(row)!=width or
            any(type(x) is not int for x in row) for row in raw):
        raise ValueError(f'{name} must be rectangular exact integer rows')
    for row in raw:
        for x in row:context.limit('max_rational_bits',abs(x).bit_length())
    return raw


def search(inputs, params, context):
    check_params(params, ('bodies','body_squared_scale','body_norm_squared',
        'heads','head_denominator','incumbent_bodies','incumbent_heads',
        'body_classes','fixed_double_body_bounds','fixed_body_norms',
        'fixed_heads','fixed_head_denominator','pure_heads','screen_epsilon',
        'max_removed','replacement_test_limit'))
    if inputs:raise ValueError('body_head_exchange takes explicit finite inputs')
    bodies=params.get('bodies');heads=params.get('heads')
    if not isinstance(bodies,list) or not bodies or not isinstance(bodies[0],list) or not bodies[0]:
        raise ValueError('Nonempty integer bodies required')
    if not isinstance(heads,list) or not heads or not isinstance(heads[0],list) or not heads[0]:
        raise ValueError('Nonempty rational-head numerator rows required')
    n,d=len(bodies),len(bodies[0]);h=len(heads[0]);palette=len(heads)
    max_removed=integer(params.get('max_removed',1),'max_removed')
    if max_removed>2:raise ValueError('Only one-for-two and two-for-three are implemented')
    test_limit=integer(params.get('replacement_test_limit',min(512,context.budget['max_trials'])),'replacement_test_limit')
    context.limit('max_dimension',d+h)
    context.limit('max_trials',max(palette,test_limit));context.limit('max_codewords',n*palette)
    D=integer(params.get('head_denominator',1),'head_denominator')
    S=integer(params.get('body_squared_scale',1),'body_squared_scale')
    R=integer(params.get('body_norm_squared'),'body_norm_squared')
    fd=integer(params.get('fixed_head_denominator',1),'fixed_head_denominator')
    for x in (D,S,R,fd):context.limit('max_rational_bits',x.bit_length())
    incumbent=params.get('incumbent_bodies')
    if not isinstance(incumbent,list) or not incumbent or any(type(x) is not int or not 0<=x<n for x in incumbent):
        raise ValueError('incumbent_bodies must be valid body indices')
    m=len(incumbent)
    bounds=params.get('fixed_double_body_bounds',[[]]);fixed=params.get('fixed_heads',[])
    norms=params.get('fixed_body_norms',[]);classes=params.get('body_classes',[0]*n)
    if not isinstance(norms,list):raise ValueError('fixed_body_norms must be a list')
    for x in norms:
        integer(x,'fixed body norm',minimum=0)
        context.limit('max_rational_bits',x.bit_length())
    f=len(norms)
    if not isinstance(classes,list) or len(classes)!=n or not isinstance(bounds,list) or not bounds or any(
            type(x) is not int or not 0<=x<len(bounds) for x in classes):
        raise ValueError('Invalid body class/envelope indices')
    context.limit('max_materialized_points',n*d+(palette+m+f)*h+len(bounds)*f)
    context.limit('max_gram_entries',n*m)
    context.limit('max_rank_work',n*m*d+palette*len(bounds)*max(f,1))
    _matrix(bodies,d,'bodies',context);_matrix(heads,h,'heads',context)
    old=params.get('incumbent_heads');_matrix(old,h,'incumbent_heads',context)
    if len(old)!=m:raise ValueError('Wrong incumbent head count')
    _matrix(fixed,h,'fixed_heads',context,empty=True)
    if len(fixed)!=f:raise ValueError('Wrong fixed head count')
    _matrix(bounds,f,'fixed_double_body_bounds',context)
    pure=params.get('pure_heads',[]);_matrix(pure,h,'pure_heads',context,empty=True)
    eps=params.get('screen_epsilon',1e-10)
    if type(eps) not in (int,float) or not isfinite(eps) or not 0<=eps<=1e-6:
        raise ValueError('screen_epsilon must be finite and in[0,1e-6]')
    if len(set(map(tuple,bodies)))!=n or len(set(map(tuple,heads)))!=palette:
        raise ValueError('Duplicate body coordinates or head palette rows')
    existing={(b,tuple(z)) for b,z in zip(incumbent,old)}
    if len(existing)!=m:raise ValueError('Duplicate incumbent body/head points')
    bound=d*max(abs(x) for row in bodies for x in row)**2
    if 2*bound>=2**62:raise ValueError('Body integer-product bound exceeds the supported proposal backend')
    X=np.array(bodies,np.int64)
    if np.any(np.sum(X*X,axis=1)!=S*R):raise ValueError('Incorrect common body norm')
    fixedN=[r*fd*fd+sum(x*x for x in z) for r,z in zip(norms,fixed)]
    pureN=[sum(x*x for x in z) for z in pure]
    if any(x<=0 for x in fixedN+pureN):raise ValueError('Zero fixed or pure point')
    oldN=[R*D*D+sum(x*x for x in z) for z in old]
    # Bounded integer sums are exactly representable in binary32. The larger
    # int64 path is exact too. The stored double body products must be integral.
    use_float=bound<2**24
    dtype=np.int8 if 2*R<=127 else (np.int16 if 2*R<=32767 else np.int64)
    dots=np.empty((n,m),dtype=dtype);Y=X[np.array(incumbent)]
    left=X.astype(np.float32) if use_float else X
    right=Y.astype(np.float32).T if use_float else Y.T
    for at in range(0,n,512):
        raw=(left[at:at+512]@right).astype(np.int64)
        if np.any((2*raw)%S):raise ValueError('Double body products are not integral')
        scaled=(2*raw)//S
        if np.any(abs(scaled)>2*R):raise ValueError('Body product exceeds the norm bound')
        dots[at:at+len(raw)]=scaled
    try:
        H=np.array(heads,float)/float(D);O=np.array(old,float)/float(D)
        F=np.array(fixed,float).reshape((f,h))/float(fd);P=np.array(pure,float).reshape((-1,h))
        beta=np.array(bounds,float);old_nf=R+np.sum(O*O,axis=1)
        fixed_nf=np.array(norms,float)+np.sum(F*F,axis=1)
        head_nf=R+np.sum(H*H,axis=1)
    except (OverflowError,ValueError) as error:
        raise ValueError('Head/envelope values are outside the floating proposal range') from error
    if any(not np.isfinite(v).all() for v in (H,O,F,P,beta,old_nf,fixed_nf,head_nf)):
        raise ValueError('Head/envelope values are outside the floating proposal range')
    classes=np.array(classes,int);groups=defaultdict(list);tested=screened=0;truncated=False
    started=monotonic()

    def accept(p,a,b):return p<=0 or p*p<=a*b

    def exact(vertex,remove):
        b,i=vertex;z=heads[i];N=R*D*D+sum(x*x for x in z)
        for beta_j,w,M in zip(bounds[int(classes[b])],fixed,fixedN):
            p=beta_j*D*fd+2*sum(x*y for x,y in zip(z,w))
            if not accept(p,N,M):return False
        for j,(w,M) in enumerate(zip(old,oldN)):
            p=int(dots[b,j])*D*D+2*sum(x*y for x,y in zip(z,w))
            if j not in remove and not accept(p,N,M):return False
        return all(accept(2*sum(x*y for x,y in zip(z,w)),N,M) for w,M in zip(pure,pureN))

    def pair(a,b):
        i,k=a;j,l=b;z,w=heads[k],heads[l]
        p=2*sum(x*y for x,y in zip(bodies[i],bodies[j]))*D*D+2*S*sum(x*y for x,y in zip(z,w))
        N=R*D*D+sum(x*x for x in z);M=R*D*D+sum(x*x for x in w)
        return p<=0 or p*p<=S*S*N*M

    def emit(vertices,removed):
        data={'backend':'finite_body_head_exchange','removed_indices':removed,
              'added_body_indices':[b for b,i in vertices],
              'added_head_indices':[i for b,i in vertices],
              'added_head_numerators':[heads[i] for b,i in vertices],
              'head_denominator':D,'input_sha256':digest(params),
              'search':{'body_count':n,'head_count':palette,'screened_candidates':screened,
                        'tested_replacement_pairs':tested,'seconds':monotonic()-started,
                        'max_removed':max_removed,'replacement_test_limit':test_limit,
                        'replacement_search_truncated':truncated,
                        'body_product_backend':'bounded exact binary32' if use_float else 'int64',
                        'body_product_absolute_bound':bound,
                        'screening':'Floating point proposals may miss feasible exchanges; every emitted local condition is checked with Python integers'},
              'scope':'Supplied incumbent geometry and fixed-family envelopes are assumptions. Independent full construction replay is required.'}
        facts={'incumbent_points':m,'added_points':len(vertices),'removed_points':len(removed),
               'gain':len(vertices)-len(removed),'new_sota':False}
        return Artifact('exchange_proposal',data,facts,Verification('unverified',
            'Exact local exchange gates under supplied family envelopes; full source geometry is not certified'))

    for i,z in enumerate(H):
        N=head_nf[i]
        if not np.all(2*(P@z)<=np.sqrt(N*np.sum(P*P,axis=1))+eps):continue
        allowed=np.all(beta+2*(F@z)[None,:]<=np.sqrt(N*fixed_nf)[None,:]+eps,axis=1)
        ids=np.flatnonzero(allowed[classes])
        threshold=np.floor(np.sqrt(N*old_nf)-2*(O@z)+eps)
        for at in range(0,len(ids),128):
            bs=ids[at:at+128];bad=dots[bs]>threshold[None,:];counts=bad.sum(axis=1)
            for b,row,count in zip(bs,bad,counts):
                b=int(b)
                if (b,tuple(heads[i])) in existing:continue
                screened+=1
                if count==0 and exact((b,i),set()):return emit([(b,i)],[])
                if 0<count<=max_removed:groups[tuple(map(int,np.flatnonzero(row)))].append((b,i))
    for removed,vertices in groups.items():
        if len(removed)!=1:continue
        for a,b in combinations(vertices,2):
            if tested>=test_limit:
                truncated=True;return emit([],[])
            tested+=1
            if pair(a,b) and exact(a,removed) and exact(b,removed):return emit([a,b],list(removed))
    if max_removed==2:
        for removed,vertices in groups.items():
            if len(removed)!=2:continue
            vertices=vertices+groups.get((removed[0],),[])+groups.get((removed[1],),[])
            context.limit('max_rank_work',n*m*d+len(vertices)**2*(d+h))
            # Pairs of candidates have exact relative body/head products. A
            # triangle gives a compatible three-point replacement family.
            cache={}
            def compatible(i,j):
                key=(min(i,j),max(i,j))
                if key not in cache:
                    context.limit('max_gram_entries',n*m+2*(len(cache)+1))
                    cache[key]=pair(vertices[i],vertices[j])
                return cache[key]
            for i,j in combinations(range(len(vertices)),2):
                if not compatible(i,j):continue
                for k in range(j+1,len(vertices)):
                    if tested>=test_limit:
                        truncated=True;return emit([],[])
                    tested+=1
                    if compatible(i,k) and compatible(j,k):
                        proposed=[vertices[t] for t in (i,j,k)]
                        if all(exact(v,removed) for v in proposed):return emit(proposed,list(removed))
    return emit([],[])


OPERATORS={'search.body_head_exchange':search}
INPUT_KINDS={'search.body_head_exchange':()}
OUTPUT_KINDS={'search.body_head_exchange':('exchange_proposal',)}
