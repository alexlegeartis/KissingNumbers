"""Exact 45D shell-and-scalar replay. python verify45.py witness.json path/to/dim44-d4-section
Uses the pinned P48n parent Gram and the universal D4 affine population identities.
Supply --moments to rerun the universal affine identity proof in the parent verifier.
"""
import argparse, base64, hashlib, itertools, json, sys, zlib
from pathlib import Path
import numpy as np

def universal_two_anchor_count():
    """11-design proof: u=<a,x>, t=<b,x>, a.a=b.b=6, a.b=3."""
    from fractions import Fraction
    from functools import lru_cache
    from math import prod
    @lru_cache(None)
    def pair(i,j):
        if (i+j)%2: return 0
        if not i+j: return 1
        if i:
            return (6*(i-1)*pair(i-2,j) if i>1 else 0) + (3*j*pair(i-1,j-1) if j else 0)
        return 6*(j-1)*pair(0,j-2)
    def moment(i,j):
        d=i+j
        if d%2: return 0
        return Fraction(52416000 * 6**(d//2)*pair(i,j),prod(range(48,48+d,2)))
    def rise(roots,sign=1):
        out=[sign]
        for z in roots:
            out=[-z*out[0]]+[out[j-1]-z*out[j] for j in range(1,len(out))]+[out[-1]]
        return out
    P=rise(range(-3,3)); Q=rise([1,2,3],-1)
    def ev(f,x):return sum(c*x**i for i,c in enumerate(f))
    pu=sum(c*moment(i,0) for i,c in enumerate(P))
    pq=sum(c*d*moment(i,j) for i,c in enumerate(P) for j,d in enumerate(Q))
    assert (pu,pq)==(26611200,9616320) and ev(P,3)==720 and ev(Q,0)==6
    assert (pu-ev(P,6)-ev(P,-6))//720==36848
    # The only exceptions to u in {-3,...,3} are +/-a.
    # With u=3, |<x,b>| and |<a-x,b>| <=3 unless x=b or a-b.
    exc=sum(ev(P,u)*ev(Q,t) for u,t in [(6,3),(-6,-3),(3,6),(3,-3)])
    assert pq-exc==1656*ev(P,3)*ev(Q,0)

def main():
    p=argparse.ArgumentParser()
    here = Path(__file__).resolve().parent
    p.add_argument("witness", type=Path, nargs="?", default=here / "witness.json")
    p.add_argument("source", type=Path, nargs="?", default=here.parent / "dim44-d4-section")
    p.add_argument("--moments",action="store_true")
    p.add_argument("--ambient",type=Path,default=here.parent/"dim41-46-p48-heads"/"partial43"/"parent_embedding.npz")
    args=p.parse_args()
    assert __debug__
    universal_two_anchor_count()
    w=json.loads(args.witness.read_text())
    parent=args.source.resolve()
    expected = {"families.npz":"267b4959b3f938c3d99ea8d5f89cad5e1d3252cd","moments.py":"7cbcda26e432fd1eb5abe23a1d28f1e5868eedec"}
    for name, sha in expected.items():
        content = (parent / name).read_bytes()
        prefix = f"blob {len(content)}\0".encode("ascii")
        assert hashlib.sha1(prefix + content).hexdigest() == sha, name
    if args.moments:
        sys.path.insert(0,str(parent))
        from moments import verify_families
        verify_families()
    content=args.ambient.read_bytes()
    assert hashlib.sha1(f"blob {len(content)}\0".encode("ascii")+content).hexdigest()=="62ab919531863b248fb5a37f121c95a8454100b5"
    with np.load(args.ambient,allow_pickle=False) as q:
        G=q["G"].astype(np.int64); B=q["BB"].astype(np.int64)
        assert np.array_equal(B@B.T,96*G)
    assert G.shape==(48,48) and np.max(abs(G))<100 and "first_two_section" in w
    assert w["ambient_lattice"]=="P48n" 
    with np.load(parent/"families.npz", allow_pickle=False) as q:
        forms=q["affine_counts"].astype(np.int64); calib=q["calibrations"].astype(np.int64)
    if "first_two_section" in w:
        shell=json.loads((args.witness.parent / w["first_two_section"]).read_text())
        if "section_zlib_b64" in shell:
            assert shell["encoding"] == "zlib-base64-int8-row-major"
            raw=zlib.decompress(base64.b64decode(shell["section_zlib_b64"],validate=True))
            assert len(raw)==1656*48
            X=np.frombuffer(raw,dtype=np.int8).astype(np.int64).reshape(1656,48)
        else:
            X=np.asarray(shell["section_points"],dtype=np.int64)
        N=np.array(shell["normals"],dtype=np.int64)
        assert X.shape==(1656,48) and N.shape==(4,48) and np.max(np.abs(X))<128 and np.max(np.abs(N))<128
    else:
        raise AssertionError("Explicit full two-anchor slice required")
    K=np.array([[6,-3,-3,-3],[-3,6,0,0],[-3,0,6,0],[-3,0,0,6]],dtype=np.int64)
    adj=np.array([[216,108,108,108],[108,108,54,54],[108,54,108,54],[108,54,54,108]],dtype=np.int64)
    assert np.array_equal(N@G@N.T,K) and np.array_equal(K@adj,324*np.eye(4,dtype=np.int64))
    assert len(X)==1656 and len(np.unique(X,axis=0))==1656
    profiles=(X@G)@N.T
    assert np.all(((X@G)*X).sum(axis=1)==6) and np.all(profiles[:,0]==3) and np.all(profiles[:,1]==0)
    counts=[int(np.all(profiles==(v if v[0]==3 else -v),axis=1).sum()) for v in calib]
    assert counts==w["calibrations"]
    # Every nearest dual representative has norm <=3 by this enumeration; if a shorter
    # representative existed it would have |v_i| <=sqrt(K_ii*3)<5, hence occurs here.
    V=np.array(list(itertools.product(range(-5,6),repeat=4)),dtype=np.int64)
    W=V@adj; sq=(W*V).sum(axis=1); labels,ix=np.unique(W%324,axis=0,return_inverse=True)
    assert len(labels)==324
    lookup={tuple(c):i for i,c in enumerate(labels)}
    opp=np.array([lookup[tuple((-c)%324)] for c in labels])
    _,mar=np.unique(np.minimum(np.arange(324),opp),return_inverse=True)
    short=np.full(324,10**8,dtype=np.int64);np.minimum.at(short,ix,sq)
    assert max(short)<=3*324 and all(short%27==0)
    param=np.array([1,*counts],dtype=np.int64); active=[]
    S=w["head_denominator"]; assert (S,w["body_scale_squared"])==(8,12)
    sizes=[]
    for c,r,h in w["families"]:
        i=lookup[tuple(c)]; q=int(r+short[i]//27)
        assert q in (72,96) and r>0
        row=int(mar[i] + (164 if q==96 else 0))
        size=int(forms[row]@param); assert size>0
        nf=r*S*S+h*h
        assert nf<=72*S*S    # two distinct bodies of one class have squared distance >= 6
        active.append((i,int(r),int(h),int(nf)));sizes.append(size)
    assert len(active)==len(set(active))
    for k,(i,r,h,n) in enumerate(active):
        for j,t,g,m in active[:k]:
            diff=lookup[tuple((labels[i]-labels[j])%324)]
            # Even nonzero projected differences lift to a nonzero parent point after adding
            # a shortest normal vector in the relevant discriminant class. For equal body
            # class and radius, coincident points must also be allowed across the two heights.
            d=0 if i==j and r==t else 72-int(short[diff]//27)
            z=(r+t-d)*S*S+2*h*g
            assert z<=0 or z*z<=n*m, (i,j,r,t,h,g,d,z*z,n*m)
    assert w["pure_head_signs"]==[-1,1] and all(4*h*h<=n for i,r,h,n in active)
    total=sum(sizes)+2
    assert total==w["total"] and total>7655213
    print(json.dumps({"dimension":45,"points":total,"old":7655213,"pure_heads":2,"gain":total-7655213,
                      "families":len(active),"cross_family_tests":len(active)*(len(active)-1)//2,
                      "D4_calibrations":counts,"parameters":param.tolist(),"sizes":sizes}))
if __name__=="__main__": main()

