"""Kerdock owner layers; exact code arithmetic, geometry and support parity."""
from itertools import combinations
from functools import reduce
from fractions import Fraction as Q
from pathlib import Path
import hashlib, importlib.util, json
HERE=Path(__file__).resolve().parent
V=HERE.parent

def pin(path,expected):
    raw=path.read_bytes()
    assert hashlib.sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()==expected, path
    return path

def load(path,name):
    p=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(p);p.loader.exec_module(m)
    return m

def span(rows):
    A=[0]
    for r in rows:A +=[w^r for w in A]
    return A

def results():
    assert __debug__
    dg=load(V/'dim65-dg-dense-owners/verify.py','dg_parent')
    P=dg.parent_module()
    raw,_=P.extended_bch(13)
    C=[sum(((a&b).bit_count()&1)<<i for i,b in enumerate(raw))for a in raw]
    D0=[reduce(int.__xor__,(r for i,r in enumerate(raw)if m>>i&1),0)for m in P.kernel(C,len(raw))]
    D,piv=P.rref(D0,64)
    assert len(D)==28 and all(x.bit_count()%4==0 for x in D)
    assert all((a&b).bit_count()%2==0 for a in D for b in D)
    # Parent BCH designed distance >=13; parity extension implies >=14;
    # divisibility by four then gives >=16 as in the pinned parent.
    X=[1]
    for _ in range(62):X.append(P.gf64_multiply(X[-1],2))
    X.append(0)
    assert len(set(X))==64
    A=[(1<<64)-1]+[sum(((x>>i)&1)<<j for j,x in enumerate(X))for i in range(6)]
    B=[A[1+i]&A[1+j]for i,j in dg.BASIS]
    assert len(P.rref(A+B,64)[1])==22 and all(P.in_span(w,D,piv)for w in A+B)
    forms={dg.form(a,0)for a in range(32)}
    assert len(forms)==32
    offsets=[0,1,2]
    groups=[{a^b for b in forms}for a in offsets]
    assert len(set.union(*groups))==96
    aff=span(A)
    assert len(aff)==128 and min(w.bit_count()for w in aff if w)==32
    diffs={a^b for a,b in combinations(forms,2)}
    assert len(diffs)==496
    dot_eval=lambda m:reduce(int.__xor__,(r for j,r in enumerate(B)if m>>j&1),0)
    mind=min((dot_eval(q)^a).bit_count() for q in diffs for a in aff)
    assert mind==28 and len(forms)*len(aff)==4096

    # Same signed axis, same axis opposite signs, and distinct axes.
    assert 27+Q(37,64)*(64-2*mind)<32
    assert 37-27<32
    assert Q(37,64)*32<32
    # Old-stock cross inequalities, with upper rational enclosures of s and r.
    s,r=Q(61,80),Q(21,4)
    assert s*s>Q(37,64) and r*r>27
    assert max(32*s,28*s+2*r)<32
    assert max(16*s,12*s+4*r)<32
    assert Q(8,3)*max(9*s,8*s+r)<32
    assert 8*s<32

    pins={66:('dim66-tail-signed-nine-refill','1f35f73a6f740051249d65e6800cdd11b65bfcca','7140a4b57643874f280267b1e297f440ce0dd0cf',32794,3200),
          67:('dim67-signed-nine-refill','e28862c4c0b20098d7b43231cb7d0c3b4db36a41','2a378b421b13115ffae15f12962bd775b2ac4d88',32669,9600)}
    ans=[]
    for d,(package,h,hc,a,b) in pins.items():
        mod=load(pin(V/package/'prep.py',h),f'stock{d}')
        W=getattr(mod,f'input{d}')()[0]
        assert len(W)==a+b and len(W)==len(set(W))
        k=[(int(x)>>64).bit_count() for x in W]
        assert k.count(0)==a and k.count(2)==b
        old=json.loads(pin(V/package/'claims.json',hc).read_text())['dimensions'][0]
        assert old['dimension']==d
        previous=old['lower_bound']; gain=(d-64)*(4096-2)
        ans.append(dict(dimension=d,previous=previous,lower_bound=previous+gain,gain=gain,
                        owner_words=(d-64)*4096,tail_poles_removed=2*(d-64)))
    assert ans==json.loads((HERE/'claims.json').read_text())['dimensions']
    return ans

if __name__=='__main__':
    print(json.dumps({'dimensions':results()},indent=2))
