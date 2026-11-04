"""Construction-program search: cyclic factor -> ring lift -> extension ->
Gray isometry -> compatible translates -> sign embedding.
Not a general theorem-discovery system. No supplied ring generator is used in synthesis.
"""
from pathlib import Path
from itertools import product
import json,importlib.util,sys
import numpy as np
import sympy as sp
P=Path(__file__).resolve().parent
REPO=P.parents[1]

def solve2(A,b):
    A=np.c_[np.array(A,dtype=np.int64)%2,np.array(b,dtype=np.int64)%2]
    r=0;piv=[]
    for j in range(A.shape[1]-1):
        ids=np.flatnonzero(A[r:,j])
        if not len(ids):continue
        k=r+ids[0];A[[r,k]]=A[[k,r]]
        for i in range(len(A)):
            if i!=r and A[i,j]:A[i]^=A[r]
        piv.append(j);r+=1
    assert not np.any((A[:,:-1].sum(axis=1)==0)&(A[:,-1]!=0))
    assert r==A.shape[1]-1, 'factor lift is not unique'
    x=np.zeros(A.shape[1]-1,dtype=int)
    for i,j in enumerate(piv):x[j]=A[i,-1]
    return x

def hensel(g,h,modulus):
    g=np.array(g,dtype=int);h=np.array(h,dtype=int)
    target=np.ones(len(g)+len(h)-1,dtype=int); m=2;history=[]
    a,b=len(g)-1,len(h)-1
    while m<modulus:
        assert np.all((target-np.convolve(g,h))%m==0)
        rhs=((target-np.convolve(g,h))//m)%2
        # Degree bounds keep both factors monic, and remove leading zero equation.
        A=np.zeros((a+b,a+b),dtype=int)
        for i in range(a):A[i:i+b+1,i]=h
        for i in range(b):A[i:i+a+1,a+i]=g
        correction=solve2(A,rhs[:-1]);g[:-1]+=m*correction[:a];h[:-1]+=m*correction[a:]
        m*=2;g%=m;h%=m
        assert np.array_equal(np.convolve(g,h)%m,target)
        history.append({'modulus':m,'g':g.tolist(),'h':h.tolist()})
    return g,h,history

def binary_factors(n):
    x=sp.Symbol('x');f=sp.Poly(sum(x**i for i in range(n)),x,modulus=2)
    _,facts=sp.factor_list(f)
    out=[]
    for g,e in facts:
        assert e==1
        h=f.exquo(g)
        if g.degree()==(n-1)//2:
            out.append(([int(g.nth(i))%2 for i in range(g.degree()+1)],
                        [int(h.nth(i))%2 for i in range(h.degree()+1)]))
    return out

def generator(g,n,q,extension=-1):
    k=n-(len(g)-1);G=np.zeros((k,n+1),dtype=np.int64)
    for i in range(k):G[i,i:i+len(g)]=g
    G[:,-1]=extension*G[:,:-1].sum(axis=1)%q
    assert np.all(np.diag(G[:,:k])%2==1) and not np.any(np.tril(G[:,:k],-1))
    return G

def gray(q):
    m=q.bit_length()-1; assert q==2**m and m>=2
    X=np.array(list(product((0,1),repeat=m-1)),dtype=int)
    table=np.array([[(a>>(m-1)&1)^sum(((a>>j)&1)*x[j] for j in range(m-1))%2 for x in X] for a in range(q)],dtype=np.uint8)
    weight=np.full(q,q//4,dtype=int);weight[0]=0;weight[q//2]=q//2
    # Generalized Gray map is not automatically a difference isometry for
    # arbitrary 2^m. Check every symbol pair before using ring weights.
    for a,b in product(range(q),repeat=2):assert np.count_nonzero(table[a]!=table[b])==weight[(a-b)%q]
    return table,weight

def enumerate_small(G,q):
    k,n=G.shape;table,weight=gray(q);hist=np.zeros(n*q//2+1,dtype=np.int64)
    for start in range(0,q**k,4096):
        labels=np.arange(start,min(start+4096,q**k),dtype=np.int64)
        coeff=(labels[:,None]//q**np.arange(k))%q
        V=coeff@G%q
        hist+=np.bincount(weight[V].sum(axis=1),minlength=len(hist))
    assert hist.sum()==q**k and hist[0]==1
    return int(np.flatnonzero(hist[1:])[0]+1),hist.tolist()

def translate_search(table,ring_length,parent_distance):
    """Search all local binary masks modulo translation stabilizer of the inner
    Gray alphabet. Lower bounds add independently across ring coordinates.
    This is conservative: it need not find all usable outer-code translates.
    """
    L=table.shape[1]
    if L>8:return {'status':'local mask search capped at 8 bits'}
    masks=np.array(list(product((0,1),repeat=L)),dtype=np.uint8)
    keys={tuple(v) for v in table};cosets=[];seen=set()
    for mask in masks:
        coset=frozenset(tuple(v^mask) for v in table)
        if coset in seen:continue
        seen.add(coset);cosets.append(mask)
    selected=[];cross=[]
    for mask in cosets:
        okay=True
        for old in selected:
            d=min(np.count_nonzero(a^mask^b^old) for a in table for b in table)
            if ring_length*d<parent_distance:okay=False;break
        if okay:selected.append(mask)
    for a in selected:
        for b in selected:
            if np.array_equal(a,b):continue
            cross.append(ring_length*min(np.count_nonzero(x^a^y^b) for x in table for y in table))
    return {'local_quotient_cosets':len(cosets),'selected_repeated_masks':[m.tolist() for m in selected],
        'multiplier_certified':len(selected),'cross_distance_lower_bound':min(cross) if cross else None,
        'selection':'greedy compatible constant-row translates; no global optimum claim'}

def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def main():
    runs=[];target_candidates=[]
    for n in (7,23):
        for index,(g0,h0) in enumerate(binary_factors(n)):
            for q in (4,8,16):
                g,h,history=hensel(g0,h0,q);G=generator(g,n,q)
                record={'cyclic_length':n,'ring_modulus':q,'factor_index':index,'ring_rank':len(G),
                    'program':['binary_cyclic_factor','coprime_Hensel_lift','negative_sum_extension','generalized_Gray','compatible_translate_union','binary_sign_embedding'],
                    'g_mod2':g0,'h_mod2':h0,'lift_history':history,'generator':G.tolist(),
                    'binary_length':(n+1)*q//2,'parent_words':q**len(G),'status':'structural_candidate_only'}
                if q**len(G)<=2**18:
                    try:table,_=gray(q)
                    except AssertionError:
                        record['status']='Gray_difference_isometry_rejected';runs.append(record);continue
                    d,hist=enumerate_small(G,q);trans=translate_search(table,n+1,d)
                    record.update(status='exhaustively_verified',minimum_distance=d,homogeneous_histogram=hist,translates=trans,
                        union_words=q**len(G)*trans['multiplier_certified'],
                        kissing_sign_layer=(4*d>=(n+1)*q//2))
                    N=(n+1)*q//2
                    if 4*d>=N:
                        geom=translate_search(table,n+1,(N+3)//4)
                        record['kissing_translate_search']=geom
                        record['kissing_sign_layer_words']=q**len(G)*geom['multiplier_certified']
                        record['kissing_distance_lower_bound']=min(d,geom['cross_distance_lower_bound'] or d)
                if n==23 and q==8:target_candidates.append(record)
                runs.append(record)
    # This comparison happens AFTER synthesis. Neither polynomial nor ring matrix
    # from the shipped certificate was input to the factor/lift algorithm.
    package=REPO/'verifications/improved/dim88-96-z8-golay'
    supplied=json.loads((package/'data/z8_generator.json').read_text())
    found=next(r for r in target_candidates if r['generator']==supplied['generator'])
    raw={'generator':found['generator'],'generator_polynomial_low_first':found['lift_history'][-1]['g'],
         'complementary_factor_low_first':found['lift_history'][-1]['h']}
    (P/'dgls_synthesized_generator.json').write_text(json.dumps(raw,indent=2))
    first=module(package/'source/replay.py','parent_replay')
    second=module(package/'source/verify_independent.py','parent_independent')
    _,a=first.verify_z8(raw);_,b=second.parent_proof(raw)
    assert a==b
    table,_=gray(8);trans=translate_search(table,24,a['binary_distance'])
    assert trans['multiplier_certified']==2
    found.update(status='two_existing_independent_exact_parent_verifiers_passed',minimum_distance=24,
                 union_words=2**37,translates=trans,kissing_sign_layer=True)
    report={'outcome':'known DGLS construction recovered within supplied construction grammar; no new code or SOTA',
            'runs':runs,'dgls_exact_receipt':a,'translation_search':trans,
            'limitations':['not autonomous invention of Hensel/Gray operations','large unverified candidates are not codes with certified distance','compatible-translate search is restricted and conservative','no new weight-enumerator recursion implemented']}
    (P/'code_search_results.json').write_text(json.dumps(report,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)))
    for r in runs:print(r['cyclic_length'],r['ring_modulus'],r['factor_index'],r['status'],r.get('minimum_distance'),r.get('union_words'),flush=True)
    print('DGLS exact verification:',a['binary_union_words'],a['binary_distance'],flush=True)
if __name__=='__main__':main()
