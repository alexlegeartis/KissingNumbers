"""The k-point moment linear program for cross-sections of an extremal even unimodular
lattice, with exact rational dual certificates.

    from kpoint_lp import run3

This is the engine behind dimensions 46, 47 (the lattice P_48) and 70, 71 (Nebe's
Gamma_72).  The mathematics is written out in ``PROOF-kpoint.md``; what follows is what
the code does.

SETTING.  L is an extremal even unimodular lattice of dimension n, minimum mu, with
N minimal vectors.  Venkov's theorem says the minimal shell, rescaled to the unit
sphere, is a spherical 11-design when n = 0 (mod 24).  Fix minimal vectors u_1..u_k with
Gram matrix G and let

    n[c_1..c_k] = #{ v in L minimal : <u_i,v> = c_i for every i }.

The cell n[0..0] counts the minimal vectors orthogonal to all of the u_i; they lie in an
(n-k)-dimensional subspace and are pairwise at angle >= 60 degrees, so

    tau(n-k) >= n[0..0].

CONSTRAINTS imposed on the unknown vector (n[c])_c:

  moments       for every exponent vector p with sum p_i <= tdes and sum p_i even,
                sum_c n[c] prod_i (c_i/mu)^{p_i} = N * (sphere moment), the sphere moment
                being evaluated exactly from the normalised Gram by Wick's theorem.
                NOTE the design property gives EVERY polynomial of degree <= 11, not only
                those with all exponents even.  Omitting the mixed ones is a real bug: it
                made two equivalent Gram matrices give different answers.  See the history
                note at the bottom of this file.
  marginals     summing out all but one coordinate must reproduce the exact one-point
                distribution M1, which is itself forced by the design property at k = 1.
  antipodal     n[c] = n[-c], because L = -L.
  translations  if <u_i,v> = mu/2 then u_i - v is again minimal, so n[c] = n[c'] with
                c'_i = mu/2 and c'_j = G_ij - c_j.  This is lattice structure that the
                moments cannot see, and it is what makes dimension 46 come out exact.
  lattice range for every lattice vector w = sum a_i u_i in the span, |v-w|^2 is a lattice
                norm, so |<v,w>| <= |w|^2/2 unless v = +-w.  Strictly stronger than the
                positive-semidefiniteness test on the Gram of (u_1..u_k,v).
  singular cells a cell whose extended Gram is singular forces v into span(u_1..u_k), so
                it has an upper bound of 1.

CERTIFICATE.  The floating-point LP is used only to FIND a dual vector y.  y is then
rounded to rationals, repaired into exact dual feasibility by lowering the coordinate of
the all-ones moment row, and the weak-duality bound

    b^T y - sum_{j : u_j finite} u_j max(0, -s_j),     s = c - A^T y,

is evaluated in exact rational arithmetic and rounded to an integer in the safe direction.
Any exactly dual-feasible y gives a valid bound, so the floating-point solver cannot
corrupt the result -- at worst it produces a weaker one.  run3 returns

    ('cert', <float LP value>, <exact certified integer>, <#cells>, <#span minima>, <#translations>)

and it is the third entry, never the second, that is claimed anywhere.

VALIDATION.  ``python kpoint_lp.py`` runs the identical code on the Leech lattice, whose
cross-sections are known independently, and on P_48:

    k=1  93150 = Lambda_23                 k=2 orthogonal   43164
    k=2  cos 1/4   44550                   k=2 60 degrees   49896 = Lambda_22
    k=3  60-degree triple  27720 = Lambda_21

all exact.

HISTORY NOTE.  The dimension 68-69 package used to carry a byte-identical copy of this
file as ``scripts/kgram9.py``; it now imports this one, so there is a single engine.
Earlier revisions imposed
only the moments in which every exponent is even -- 56 equations for k = 3 instead of 161.
The bug was visible in the output and went unnoticed: the Gram matrices (3,-3,0) and 3*A_3
are equivalent (U^T G U = 3*A_3 for an integer U of determinant -1) yet gave 6462480 and
5194969.  With the mixed moments both give 6687320.  Any future version of this machinery
should be checked on a pair of equivalent Gram matrices; it is cheap and decisive.
"""

from fractions import Fraction as F
from math import gcd
import os
import numpy as np, itertools, sys, time
from scipy.optimize import linprog
def matchings(items):
    if not items:
        yield []
        return
    a = items[0]
    for i in range(1, len(items)):
        rest = items[1:i] + items[i + 1:]
        for r in matchings(rest):
            yield [(a, items[i])] + r


def wick(mult, Gn):
    """Gaussian moment E[prod <a_i,Z>^{mult_i}] as a sum over perfect matchings; exact."""
    items = []
    for i, m in enumerate(mult):
        items += [i] * m
    if len(items) % 2:
        return F(0)
    tot = F(0)
    for mt in matchings(items):
        p = F(1)
        for a, b in mt:
            p *= Gn[a][b]
        tot += p
    return tot

def latrange(G,mu,K,AMAX=3,NUMAX=48):
    """lattice vectors w = sum a_i u_i in the span, with their norms"""
    out=[]
    for a in itertools.product(range(-AMAX,AMAX+1),repeat=K):
        if all(x==0 for x in a): continue
        nu=sum(a[i]*a[j]*G[i][j] for i in range(K) for j in range(K))
        if 0<nu<=NUMAX: out.append((a,nu))
    return out

def allowed(cs,W,mu):
    """is the cell cs compatible with every lattice vector in the span?"""
    for a,nu in W:
        c=sum(a[i]*cs[i] for i in range(len(cs)))
        if abs(c)*2>nu:
            if not (nu==mu and abs(c)==mu): return False
    return True


def _integralise(rows,rhs):
    """Scale every constraint to integers and divide out the gcd.

    Row i has entries (c/mu)^p and a right-hand side carrying the sphere-moment denominator,
    so the two share a denominator L_i.  Replacing (A_i,b_i) by (L_i A_i, L_i b_i) leaves the
    feasible set unchanged and rescales the dual as y_i -> y_i/L_i, so no bound is affected.
    Rows are returned sparse: only the moment rows are dense, the rest have a few nonzeros."""
    A=[];b=[];L=[]
    for r,rh in zip(rows,rhs):
        den=1
        for x in r: den=den*x.denominator//gcd(den,x.denominator)
        den=den*rh.denominator//gcd(den,rh.denominator)
        ent=[(k,int(x*den)) for k,x in enumerate(r) if x]
        bi=int(rh*den); g=abs(bi)
        for _,v in ent: g=gcd(g,abs(v))
        if g>1:
            ent=[(k,v//g) for k,v in ent]; bi//=g; den=F(den,g)
        A.append(ent); b.append(bi); L.append(F(den))
    return A,b,L

def _cert_integer(rows,rhs,ub,m,idx,tags,Z,sense,y):
    """Exact integer lower bound on x[Z] from the numeric dual y.

    For  min c.x  s.t.  Ax=b, 0<=x<=u,  ANY vector y gives the valid bound

        y.b - sum_j u_j max(0,-slack_j),      slack_j = c_j - (A^T y)_j,

    as long as slack_j >= 0 wherever x_j is unbounded (weak duality; the bounded j are paid
    for by the penalty).  With A,b integral and y carried as a/D over ONE common denominator,

        D*slack_j = D*c_j - sum_i a_i A_ij

    is an integer, so feasibility is an integer comparison and the bound is a single division.
    The previous version instead carried y as unrelated Fractions from limit_denominator(D):
    against row denominators up to mu^tdes the running denominator of one inner product grew
    past 10^28, so every slack cost a gcd on 100-bit integers and the check took over an hour.
    Nothing in the mathematics needed it -- the certificate here is the same inequality."""
    A,b,L=_integralise(rows,rhs)
    nr=len(A); SG=int(sense); zj=idx[Z]
    free=[k for k in range(m) if ub[k] is None]
    bdd=[(k,int(ub[k])) for k in range(m) if ub[k] is not None]
    ai=next(i for i,t in enumerate(tags) if t[0]=='mom' and all(q==0 for q in t[1]))
    assert len(A[ai])==m, 'the p=0 moment row must have a nonzero on every cell'
    e0=A[ai][0][1]                      # all its entries are equal, so shifting y[ai] moves
    assert e0>0                         # every slack by the same amount
    best=None
    for D in (10**6,10**9,10**12,10**15,10**18):
        a=[int(round(float(y[i])/float(L[i])*D)) for i in range(nr)]
        S=[0]*m
        for i in range(nr):
            ci=a[i]
            if ci:
                for k,v in A[i]: S[k]+=ci*v
        T=[(D*SG if k==zj else 0)-S[k] for k in range(m)]
        d=min(T[k] for k in free)
        if d<0:                         # Lemma 3 repair: sh = floor(d/e0), so sh*e0 <= d
            sh=d//e0; a[ai]+=sh; se=sh*e0
            for k in range(m): T[k]-=se
            d=min(T[k] for k in free)
        if d<0: continue
        num=sum(a[i]*b[i] for i in range(nr))
        for k,u in bdd:
            if T[k]<0: num-=u*(-T[k])
        # ceil for a minimum; for a maximum the tight integer bound is floor(-B), not
        # -floor(B) -- the two differ by one whenever B is not an integer
        v=-((-num)//D) if sense>0 else (-num)//D
        if best is None or (v>best if sense>0 else v<best): best=v
    return best

def _cert_rational(rows,rhs,ub,m,idx,tags,Z,sense,y):
    """The original rational certificate, kept for comparison (see _cert_integer)."""
    SG=F(sense)
    def slack(yv): return [(SG if j==idx[Z] else F(0))-sum(yv[i]*rows[i][j] for i in range(len(rows))) for j in range(m)]
    free=[j for j in range(m) if ub[j] is None]
    ai=next(i for i,t in enumerate(tags) if t[0]=='mom' and all(p==0 for p in t[1]))
    best=None
    for DEN in (10**6,10**9,10**12,10**15,10**18):
        yq=[F(float(v)).limit_denominator(DEN) for v in y]
        sl=slack(yq); d=min(sl[j] for j in free)
        if d<0: yq[ai]+=d; sl=slack(yq)
        if min(sl[j] for j in free)<0: continue
        pen=sum(F(int(ub[j]))*max(F(0),-sl[j]) for j in range(m) if ub[j] is not None)
        bnd=sum(yq[i]*rhs[i] for i in range(len(rows)))-pen
        v=-((-bnd.numerator)//bnd.denominator)
        if sense<0: v=(-bnd.numerator)//bnd.denominator      # floor(-B), see _cert_integer
        if best is None or (v>best if sense>0 else v<best): best=v
    return best

def build_system(n,N,mu,IPS,M1,G,AMAX=2,tdes=11):
    """The k-point system, exactly as run3 has always built it.

    Returned so that a certifier other than run3's own can work on the SAME rows:
    a bound that rests on a rebuilt system rests on two implementations agreeing,
    which is a weaker thing to check than one implementation being right.
    """
    K=len(G)
    Gn=[[F(G[i][j],mu) for j in range(K)] for i in range(K)]
    def sphmom(ps):
        M=sum(ps); den=F(1)
        for j in range(M//2): den*=(n+2*j)
        return wick(list(ps),Gn)/den
    vals=[0]+[s*c for c in IPS for s in (1,-1)]+[mu,-mu]
    def det(cs):
        Mx=[[F(G[i][j]) for j in range(K)]+[F(cs[i])] for i in range(K)]
        Mx.append([F(c) for c in cs]+[F(mu)])
        s=K+1; d=F(1); Mx=[r[:] for r in Mx]
        for i in range(s):
            p=None
            for r in range(i,s):
                if Mx[r][i]!=0: p=r; break
            if p is None: return F(0)
            if p!=i: Mx[i],Mx[p]=Mx[p],Mx[i]; d=-d
            d*=Mx[i][i]; inv=Mx[i][i]
            for r in range(i+1,s):
                f=Mx[r][i]/inv
                if f: Mx[r]=[a-f*b for a,b in zip(Mx[r],Mx[i])]
        return d
    # minimal vectors in the span: w = sum a_i u_i with |w|^2 = mu
    spanmin=[]
    for a in itertools.product(range(-AMAX,AMAX+1),repeat=K):
        if all(x==0 for x in a): continue
        nrm=sum(a[i]*a[j]*G[i][j] for i in range(K) for j in range(K))
        if nrm!=mu: continue
        prof=tuple(sum(a[i]*G[i][l] for i in range(K)) for l in range(K))
        spanmin.append((a,prof))
    cells=[];known={};rig=[]
    profs={p for _,p in spanmin}
    W=latrange(G,mu,K)
    for cs in itertools.product(vals,repeat=K):
        d=det(cs)
        if d<0: continue
        if not allowed(cs,W,mu): continue
        if cs in profs: known[cs]=1
        else:
            cells.append(cs)
            if d==0: rig.append(cs)
    idx={c:i for i,c in enumerate(cells)}; m=len(cells)
    rows=[];rhs=[];tags=[]
    # Venkov: the shell is a tdes-design, so EVERY polynomial of degree <= tdes integrates
    # exactly -- not merely those with all exponents even.  Iterate over all exponent vectors
    # with even total degree <= tdes; odd total degree vanishes by antipodal symmetry.  For
    # k = 3 this is 161 equations rather than 56, and it is what pins Ozeki's tau in [184,278].
    for ps in itertools.product(range(tdes+1),repeat=K):
        t=sum(ps)
        if t>tdes or t%2: continue
        if t==0 and any(ps): continue
        r=[]
        for cs in cells:
            v=F(1)
            for c,p in zip(cs,ps):
                if p: v*=F(c,mu)**p
            r.append(v)
        ex=F(0)
        for cs,cnt in known.items():
            v=F(1)
            for c,p in zip(cs,ps):
                if p: v*=F(c,mu)**p
            ex+=cnt*v
        rows.append(r); rhs.append(F(N)*sphmom(ps)-ex); tags.append(('mom',ps))
    for i in range(K):
        for a in vals:
            if abs(a)==mu: continue
            r=[F(1) if c[i]==a else F(0) for c in cells]
            ex=sum(cnt for cs,cnt in known.items() if cs[i]==a)
            rows.append(r); rhs.append(F(M1[abs(a)])-ex); tags.append(('marg',i,a))
    for c in cells:
        d=tuple(-x for x in c)
        if c<d and d in idx:
            r=[F(0)]*m; r[idx[c]]=F(1); r[idx[d]]=F(-1); rows.append(r); rhs.append(F(0)); tags.append(('sym',c))
    h=mu//2; nt=0
    for a,prof in spanmin:
        for c in cells:
            ip=sum(a[i]*c[i] for i in range(K))
            if ip!=h: continue
            cp=tuple(prof[l]-c[l] for l in range(K))
            if cp==c: continue
            if cp in idx:
                r=[F(0)]*m; r[idx[c]]=F(1); r[idx[cp]]=F(-1)
                rows.append(r); rhs.append(F(0))
            elif cp in known:
                r=[F(0)]*m; r[idx[c]]=F(1); rows.append(r); rhs.append(F(known[cp]))
            else:
                r=[F(0)]*m; r[idx[c]]=F(1); rows.append(r); rhs.append(F(0))
            tags.append(('tr',a,c)); nt+=1
    ub=[None]*m
    for c in rig: ub[idx[c]]=1.0
    return dict(K=K,cells=cells,idx=idx,rows=rows,rhs=rhs,ub=ub,tags=tags,
                known=known,spanmin=spanmin,nt=nt,m=m,rig=rig,vals=vals)
def run3(n,N,mu,IPS,M1,G,target=None,AMAX=2,verbose=False,tdes=11,sense=1,
         extra_solves=False):
    _S=build_system(n,N,mu,IPS,M1,G,AMAX,tdes)
    K=_S['K']; cells=_S['cells']; idx=_S['idx']; rows=_S['rows']; rhs=_S['rhs']
    ub=_S['ub']; tags=_S['tags']; known=_S['known']; spanmin=_S['spanmin']
    nt=_S['nt']; m=_S['m']
    # ---- numeric duals ------------------------------------------------------------------
    # Solve in shell fractions z = x/N so the system is well conditioned: min c.z subject to
    # Az = b/N, 0 <= z <= ub/N has the SAME dual feasible set A^T y <= c as the original, and
    # y.b is the original optimum.  The exact certificate below re-derives everything from
    # rows/rhs/ub, so nothing numeric decides anything.
    A0=np.array([[float(x) for x in r] for r in rows]); b0=np.array([float(x) for x in rhs])
    rs=np.maximum(np.abs(A0).max(axis=1),1e-300); A=A0/rs[:,None]; b=b0/rs/N
    ubn=[(None if u is None else u/N) for u in ub]
    bnds=list(zip([0.0]*m,ubn))
    Z=tuple([0]*K) if target is None else tuple(target)
    if Z not in idx: return None
    obj=np.zeros(m); obj[idx[Z]]=float(sense)   # sense=-1 maximises the cell
    budget=float(os.environ.get('KN_BUDGET','900'))
    t0=time.time()
    # A per-SOLVE limit as well as an overall one.  Checking the budget between solves is
    # not enough: the A_4 four-point program at tdes >= 9 is degenerate enough that a single
    # HiGHS call runs for tens of minutes, which is how this used to hang with no way out.
    def opts(presolve=True):
        return {'time_limit': max(5.0, budget-(time.time()-t0)), 'presolve': presolve}
    cands=[]
    # extra_solves=True adds the presolve=False variants and keeps EVERY dual instead of
    # stopping at the first success.  On these systems HiGHS's presolve returns points that
    # violate the constraints by tens of vectors, and a dual read off one of those certifies
    # up to 0.3% low -- which is the whole of the apparent spread across isometric Grams.
    # The certificate is the max over candidates and each is verified exactly, so this can
    # only raise a bound.  Off by default so that no shipped number moves.
    tries=[(m,p) for m in ('highs','highs-ipm','highs-ds')
           for p in ((False,True) if extra_solves else (True,))]
    for meth,pre in tries:
        if time.time()-t0>budget: break
        r0=linprog(c=obj,A_eq=A,b_eq=b,bounds=bnds,method=meth,options=opts(pre))
        if r0.success and getattr(r0,'eqlin',None) is not None:
            cands.append(('eq:%s/pre=%s'%(meth,pre),r0,np.asarray(r0.eqlin.marginals)/rs))
            if not extra_solves: break
    # Always try the relaxed form too, budget permitting: the certificate is the MAX over
    # candidates, so an extra dual can only raise the bound.
    if time.time()-t0<=budget:
        Au=np.vstack([A,-A]); done=False
        for EPS in (1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5):
            if done or time.time()-t0>budget: break
            bu=np.concatenate([b+EPS,-b+EPS])
            for meth in ('highs','highs-ipm','highs-ds'):
                if time.time()-t0>budget: break
                r0=linprog(c=obj,A_ub=Au,b_ub=bu,bounds=bnds,method=meth,options=opts())
                if r0.success and getattr(r0,'ineqlin',None) is not None:
                    lam=np.asarray(r0.ineqlin.marginals); hh=len(lam)//2
                    cands.append(('ub%.0e:%s'%(EPS,meth),r0,(lam[:hh]-lam[hh:])/rs))
                    done=True; break
    if not cands: return None
    if verbose:
        print('   duals: %s   (%.0fs)'%(', '.join(c[0] for c in cands),time.time()-t0))
    # ---- the exact certificate, best over the candidates --------------------------------
    cert=_cert_rational if os.environ.get('KN_CERT')=='rational' else _cert_integer
    best=None; res=cands[0][1]
    for lbl,r0,y in cands:
        v=cert(rows,rhs,ub,m,idx,tags,Z,sense,y)
        if v is None: continue
        if best is None or (v>best if sense>0 else v<best): best=v; res=r0
    if best is None: return ('lp',res.fun*N,m,len(spanmin),nt)
    return ('cert',res.fun*N,best,m,len(spanmin),nt)
if __name__ == '__main__':
    import time
    print(__doc__)
    ok = True

    def check(label, got, want):
        global ok
        good = (got == want)
        ok = ok and good
        print("   %-26s %-12s expected %-12s %s"
              % (label, got, want, "MATCH" if good else "*** MISMATCH ***"))

    print("=" * 78)
    print("VALIDATION 1: the Leech lattice, whose cross-sections are known independently")
    print("=" * 78)
    M1L = {0: 93150, 1: 47104, 2: 4600}
    for G, name, truth in (
            ([[4, 0], [0, 4]], 'k=2 orthogonal', 43164),
            ([[4, 1], [1, 4]], 'k=2 cos 1/4', 44550),
            ([[4, 2], [2, 4]], 'k=2 60 degrees', 49896),
            ([[4, 0, 0], [0, 4, 0], [0, 0, 4]], 'k=3 orthogonal', 19530),
            ([[4, 2, 2], [2, 4, 2], [2, 2, 4]], 'k=3 60-degree triple', 27720),
            ([[4, 2, 2, 2], [2, 4, 2, 2], [2, 2, 4, 2], [2, 2, 2, 4]], 'k=4 A_4', 15540),
            ([[4, -2, -2, -2], [-2, 4, 0, 0], [-2, 0, 4, 0], [-2, 0, 0, 4]], 'k=4 D_4', 17400)):
        t = time.time()
        r = run3(24, 196560, 4, [1, 2], M1L, G)
        check(name, r[2], truth)
    print("   (49896 is Lambda_22, 27720 is Lambda_21, 17400 is Lambda_20.  Every target here")
    print("    is COUNTED, not quoted: verifications/improved/dim68-69-gamma72-cross-sections/")
    print("    scripts/leech_truth.py builds all 196560 minimal vectors from the Golay code and")
    print("    counts them for twelve independent tuples of each Gram, and reproduces the")
    print("    marginals 93150/47104/4600 that this program takes as input.  The programme is")
    print("    The cross-section is NOT a function of the Gram: two orthogonal triples")
    print("    with Gram mu*I_3 give 19530 and 19962.  So the programme BRACKETS every")
    print("    embedding rather than pinning one, and on the Leech it is TIGHT -- the")
    print("    certified minimum is attained for every Gram here, and for mu*I_3 both")
    print("    ends [19530, 19962] are attained.  Gamma_72 quotes the LOWER end, which")
    print("    is a bound valid for every embedding, which is what a claim needs.)")

    print()
    print("=" * 78)
    print("VALIDATION 2: P_48, where the k = 1 answer is forced by the design property")
    print("=" * 78)
    M1P = {0: 23766960, 1: 12608784, 2: 1678887, 3: 36848}
    r = run3(48, 52416000, 6, [1, 2, 3], M1P, [[6, 3], [3, 6]])
    check('k=2 60 deg -> dim 46', r[2], 12309600)
    r = run3(48, 52416000, 6, [1, 2, 3], M1P, [[6, 3], [3, 6]], sense=-1)
    check('the same, MAXIMISED', r[2], 12309600)
    print("   (minimum = maximum, so that cross-section has exactly 12309600 minimal")
    print("    vectors; independently confirmed by Ozeki's degree-3 Siegel theta series.)")

    print()
    print("=" * 78)
    print("VALIDATION 3: Gamma_72 -> dimension 70")
    print("=" * 78)
    M1G = {0: 2603658750, 1: 1512243200, 2: 280928256, 3: 13959168, 4: 127800}
    r = run3(72, 6218175600, 8, [1, 2, 3, 4], M1G, [[8, 4], [4, 8]])
    check('k=2 60 deg -> dim 70', r[2], 1249778250)

    print()
    print("ALL CHECKS PASSED" if ok else "*** FAILURES ***")
    raise SystemExit(0 if ok else 1)
