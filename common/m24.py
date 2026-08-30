"""Automorphisms of a Golay code given by its 759 octads: backtracking using S(5,8,24)."""
import numpy as np, sys, os, random
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
def octad_bitmasks(G):
    return [int(''.join(str(int(b)) for b in w[::-1]),2) for w in G if w.sum()==8]
def build_tables(octs):
    oset=set(octs)
    five={}
    for o in octs:
        pts=[i for i in range(24) if (o>>i)&1]
        from itertools import combinations
        for c in combinations(pts,5):
            m=0
            for x in c: m|=1<<x
            five[m]=o
    byp=[[] for _ in range(24)]
    for o in octs:
        for i in range(24):
            if (o>>i)&1: byp[i].append(o)
    return oset,five,byp
def find_autos(G, n_want=8, seed=0, tries=400):
    octs=octad_bitmasks(G); oset,five,byp=build_tables(octs)
    rng=random.Random(seed)
    res=[]
    pts_of={o:[i for i in range(24) if (o>>i)&1] for o in octs}
    def ok(pi, assigned_mask, t):
        # check every octad containing t that has >=5 assigned points
        for o in byp[t]:
            am=o & assigned_mask
            k=bin(am).count('1')
            if k<5: continue
            idx=[i for i in range(24) if (am>>i)&1]
            m5=0
            for x in idx[:5]: m5|=1<<pi[x]
            tgt=five.get(m5)
            if tgt is None: return False
            for x in idx:
                if not ((tgt>>pi[x])&1): return False
        return True
    def bt(pi, used, t):
        if t==24: return True
        order=list(range(24)); rng.shuffle(order)
        for v in order:
            if used>>v & 1: continue
            pi[t]=v
            if ok(pi, (1<<(t+1))-1, t) and bt(pi, used|(1<<v), t+1): return True
            pi[t]=-1
        return False
    for _ in range(tries):
        pi=[-1]*24
        if bt(pi,0,0):
            p=np.array(pi,dtype=np.int64)
            if not any((p==q).all() for q in res): res.append(p)
            if len(res)>=n_want: break
    return res
if __name__ == "__main__":
    # Self-test: find automorphisms of the Golay code (elements of M_24) and check that
    # each really permutes the code.  The Golay code is built here from golay.py rather
    # than loaded from anywhere, so this runs standalone.
    import time
    from golay import golay24
    G = golay24()
    t0 = time.time()
    A = find_autos(G, n_want=int(sys.argv[1]) if len(sys.argv) > 1 else 6, seed=1)
    print("found %d automorphisms in %.0fs" % (len(A), time.time() - t0))
    GS = set(map(tuple, G.tolist()))
    ok = True
    for p in A[:6]:
        good = all(tuple(w[p].tolist()) in GS for w in G[[1, 2, 4, 8, 16, 32]])
        ok = ok and good
        print("   perm ok: %s  %s ..." % (good, p[:8].tolist()))
    print("ALL CHECKS PASSED" if ok else "*** FAILURES ***")
    raise SystemExit(0 if ok else 1)
