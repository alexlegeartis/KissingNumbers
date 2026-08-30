"""Gamma_72 on the GPU: pool generation and disjoint-class building, in EXACT arithmetic.

WHAT THIS IS FOR.  Dimensions 87-95 use the cap construction over Gamma_72, whose gain is
    4 * (sum of the |T| largest sizes among pairwise-DISJOINT classes) + 2 * (next |P|),
a class being a set of minimal LINES with pairwise |<u,v>| <= 2.  So the whole question is
how many minimal lines can be packed into ~31 050 disjoint classes.  Greedy on a candidate
set of n lines returns about ln(1 + n p)/p lines, where p is the density of conflicting
pairs, so BOTH n AND p matter and nothing else does:

    one automorphism orbit   n = 2 057 093   p = 0.010241   ->   972
    orbit + difference rule  n = 3 387 346   p = 0.005469   ->  1797   (record 2106 w/ search)
    this module, target      n = 80 000 000  p ~ 0.0055     ->  2375

and 31 050 classes drawn that way total 69.7M lines against the present family's 46.2M.

EXACTNESS.  <x,y> = (x G) . y with integer x, y, G, and the two ranges are BOUNDS, not
measurements.  Gamma_72 is unimodular, so G^-1 is integral, and for a minimal x

    |x_j| = |<x, b*_j>| <= |x| |b*_j| = sqrt(8 * (G^-1)_jj) = sqrt(8 * 72) = 24 ,
    |(x G)_j| = |<x, b_j>| <= |x| |b_j| = sqrt(8 * G_jj)    = sqrt(8 *  8) =  8 ,

max_j (G^-1)_jj = 72 and max_j G_jj = 8 being read off the Gram.  Hence every 72-term dot
product is an integer of modulus at most 72 * 24 * 8 = 13 824.  (An earlier version assumed
16 and 4 from a sample of the working pool; the walk immediately produced a coordinate of 20.)
Three back ends are tried in order and each is checked against an int64 reference at start-up:

  * torch._int_mm  -- int8 x int8 -> int32 on the T4's INT8 tensor cores, exact by definition;
  * float16 in, float32 accumulate -- the inputs (<= 24 and <= 8) are exact in float16, the
    accumulation is exact in float32 (13 824 << 2^24), and only the store back to float16
    rounds: values below 2048 are exact and a value of modulus >= 2048 comes back to within 8
    of the truth, so it can never be mistaken for one of modulus <= 2, the only test made;
  * float32 -- exact to 2^24, a factor 3600 of headroom.

TF32 MUST BE OFF for the float32 path: it keeps 10 mantissa bits (exact only to 1024) and
would silently corrupt inner products of size up to 4608.  Both switches are off below.

DETERMINISM, AND ITS LIMIT.  Every random CHOICE is made by numpy's PCG64 from the run's seed
and only then moved to the device; the GPU does nothing but exact integer GEMM.  That is what
lets the certificate be 100 kB instead of the 5 GB the lines themselves would take -- but it is
not the same as reproducibility on any machine.  The pool is ASSEMBLED by torch primitives whose
tie-breaking is unspecified (dedup_into sorts by a 64-bit key), and its ROW ORDER is what the
greedy window is drawn from, so two software stacks can produce two different pools from the
same seed.  Measured: a certificate made on a T4 under one Colab image did not reproduce under a
later one, while two builds inside one process agreed bit for bit.  Everything downstream is
indexed by pool POSITION, so pool_digest() records a fingerprint, the state carries it, and a
--resume against a pool that does not match is refused rather than silently reusing lines.
"""
import numpy as np, torch, os, time, json, hashlib
import fused as _fused

torch.backends.cuda.matmul.allow_tf32 = False
torch.backends.cudnn.allow_tf32 = False
torch.backends.cuda.matmul.allow_fp16_reduced_precision_reduction = False

HERE = os.path.dirname(os.path.abspath(__file__))
EXACT_LIMIT = 2 ** 24
MAXC = 24          # sqrt(8 * max_j (G^-1)_jj) = sqrt(8*72), exact bound on |x_j|
MAXG = 8           # sqrt(8 * max_j G_jj)      = sqrt(8* 8), exact bound on |(x G)_j|

_BACKEND = None


def device():
    return torch.device('cuda' if torch.cuda.is_available() else 'cpu')


def load_gram():
    return np.load(os.path.join(HERE, 'gamma72_gram.npy')).astype(np.int64)


def check_bounds():
    """recompute MAXC and MAXG from the Gram itself, so the constants can never drift"""
    G = load_gram()
    Gi = np.rint(np.linalg.inv(G.astype(np.float64))).astype(np.int64)
    assert (Gi @ G == np.eye(72, dtype=np.int64)).all(), "Gamma_72 Gram is not unimodular"
    c = int(np.floor(np.sqrt(8.0 * int(Gi.diagonal().max()))))
    g = int(np.floor(np.sqrt(8.0 * int(G.diagonal().max()))))
    assert c == MAXC and g == MAXG, "coordinate bounds are %d, %d, not %d, %d" % (c, g, MAXC, MAXG)
    assert 72 * MAXC * MAXG < EXACT_LIMIT
    return 72 * MAXC * MAXG


def load_gens():
    return np.load(os.path.join(HERE, 'gamma72_gens.npy')).astype(np.int64)


def canon_t(V):
    """sign-normalise rows of an int8 GPU tensor so the first non-zero coordinate is positive"""
    nz = V != 0
    first = torch.argmax(nz.to(torch.uint8), dim=1)
    lead = V.gather(1, first[:, None]).squeeze(1)
    s = torch.where(lead < 0, torch.tensor(-1, dtype=V.dtype, device=V.device),
                    torch.tensor(1, dtype=V.dtype, device=V.device))
    return V * s[:, None]


def canon(V):
    nz = V != 0
    first = np.argmax(nz, axis=1)
    s = np.sign(V[np.arange(len(V)), first]).astype(V.dtype)
    return V * s[:, None]


# ---------------------------------------------------------------- exact GEMM back ends

def _mm_int8(A, B):
    return torch._int_mm(A, B)


def _mm_f16(A, B):
    return (A.to(torch.float16) @ B.to(torch.float16)).to(torch.int32)


def _mm_f32(A, B):
    return (A.to(torch.float32) @ B.to(torch.float32)).to(torch.int32)


def select_backend(dev, verbose=True):
    """pick the fastest back end that reproduces an int64 reference exactly, and say which"""
    global _BACKEND
    b = check_bounds()
    if verbose:
        print("   |<x,y>| <= 72 * %d * %d = %d < 2^24, from the Gram" % (MAXC, MAXG, b), flush=True)
    rng = np.random.default_rng(0)
    A = rng.integers(-MAXC, MAXC + 1, size=(256, 72)).astype(np.int8)
    B = rng.integers(-MAXG, MAXG + 1, size=(72, 256)).astype(np.int8)
    ref = A.astype(np.int64) @ B.astype(np.int64)
    At, Bt = torch.tensor(A, device=dev), torch.tensor(B, device=dev)
    for name, fn in (('int8', _mm_int8), ('fp16', _mm_f16), ('fp32', _mm_f32)):
        try:
            got = fn(At, Bt).cpu().numpy().astype(np.int64)
        except Exception as e:                                   # unsupported on this device
            if verbose:
                print("   backend %s unavailable (%s)" % (name, type(e).__name__), flush=True)
            continue
        if (got == ref).all():
            _BACKEND = fn
            if verbose:
                print("   GEMM backend: %s (verified against int64)" % name, flush=True)
            return name
        if verbose:
            print("   backend %s INEXACT -- rejected" % name, flush=True)
    raise RuntimeError("no exact GEMM backend on this device")


def assert_fp16_exact(dev):
    """the float16-in / float32-accumulate path used by build_classes reproduces int64 on the
    only test it makes, |ip| <= 2, for every attainable inner product"""
    rng = np.random.default_rng(1)
    A = rng.integers(-MAXC, MAXC + 1, size=(1024, 72))
    B = rng.integers(-MAXG, MAXG + 1, size=(72, 1024))
    ref = A @ B
    got = (torch.tensor(A, dtype=torch.float16, device=dev)
           @ torch.tensor(B, dtype=torch.float16, device=dev)).to(torch.float64).cpu().numpy()
    assert ((np.abs(ref) <= 2) == (np.abs(got) <= 2)).all(), "float16 path misclassifies |ip|<=2"
    assert (got[np.abs(ref) < 2048] == ref[np.abs(ref) < 2048]).all(), "float16 path inexact"


def mm(A, B):
    """exact integer GEMM, with the int8 back end's shape rules handled here.

    torch._int_mm insists that the inner and output dimensions be multiples of 8 and the row
    count be at least 32 (on CUDA; the CPU implementation does not enforce it, which is why
    this only appeared on the T4 -- the accepted block of a greedy round has whatever width
    the block greedy returned, 73 on the run that caught it).  Padding with zeros cannot
    change the product of the real entries, so pad and slice.
    """
    if _BACKEND is None:
        select_backend(A.device, verbose=False)
    if _BACKEND is not _mm_int8:
        return _BACKEND(A, B)
    m, k = A.shape
    n = B.shape[1]
    pm, pk, pn = (-m) % 32, (-k) % 8, (-n) % 8
    if pm or pk or pn:
        A = torch.nn.functional.pad(A, (0, pk, 0, pm))
        B = torch.nn.functional.pad(B, (0, pn, 0, pk))
    return _mm_int8(A.contiguous(), B.contiguous())[:m, :n]


# ---------------------------------------------------------------- hashing / dedup

def _hash_keys(X, r1, r2):
    """one int64 key per row, exactly equal for equal rows; distinct rows collide with
    probability ~n^2/2^65, and a collision can only DISCARD a line, never keep a duplicate"""
    n = len(X)
    h1 = torch.zeros(n, dtype=torch.int64, device=X.device)
    h2 = torch.zeros(n, dtype=torch.int64, device=X.device)
    for j in range(72):
        c = X[:, j].to(torch.int64)
        h1 += c * int(r1[j])
        h2 += c * int(r2[j])
    return (h1 & 0xFFFFFFFF) * (1 << 32) + (h2 & 0xFFFFFFFF)


def dedup_into(X, keys, Xout, Hs, nout, verbose=False):
    """Append the rows of X whose key is not already among the sorted keys Hs.

    Hs is kept sorted by MERGING rather than by re-sorting the concatenation.  torch.sort
    returns indices as well as values, so sorting a 95M-key array costs about 2.3 GB of
    transient allocation on top of the 6.4 GB pool; at that size the caching allocator starts
    evicting and the pool phase went from 133 s to over 15 minutes.  Merging two already
    sorted arrays is three linear passes and one boolean mask.
    """
    o = torch.argsort(keys)
    keys, X = keys[o], X[o]
    del o
    fresh = torch.ones(len(keys), dtype=torch.bool, device=X.device)
    fresh[1:] = keys[1:] != keys[:-1]                       # duplicates inside the batch
    keys, X = keys[fresh], X[fresh]
    del fresh
    if nout:
        pos = torch.searchsorted(Hs, keys).clamp_(max=nout - 1)
        seen = Hs[pos] == keys
        keys, X = keys[~seen], X[~seen]
        del pos, seen
    take = min(len(keys), len(Xout) - nout)
    if take <= 0:
        return nout, Hs
    keys, X = keys[:take], X[:take]
    Xout[nout:nout + take] = X
    if not nout:
        return take, keys
    ins = torch.searchsorted(Hs, keys) + torch.arange(take, device=keys.device)
    merged = torch.empty(nout + take, dtype=torch.int64, device=keys.device)
    merged[ins] = keys
    mask = torch.ones(nout + take, dtype=torch.bool, device=keys.device)
    mask[ins] = False
    merged[mask] = Hs
    return nout + take, merged


def residency(pool, scan, block, pchunk):
    """GPU bytes the class phase holds: the pool (int8), the float16 window and the copy the
    compaction makes of it, the alive index, and the product tile.  Returned in GiB as
    (total, pool, window, tile)."""
    G = float(2 ** 30)
    p = pool * 72 / G
    w = scan * 72 * 2 * 2 / G + pool * 4 / G
    t = pchunk * block * 2 / G
    return (p + w + t, p, w, t)


def autosize(block=128, pchunk=2097152, frac=0.60, ratio=10, cap=400_000_000):
    """Choose pool and scan for the free memory of this device.

    Keep the total under `frac` of what is free: measured on a T4, a class takes 0.23 s at
    9.2 GB of residency and 6 s at 11.5 GB out of 14.7 GB, with nothing about the arithmetic
    changed -- past roughly 70%% the caching allocator starts evicting and re-serving the same
    buffers.  scan = pool/ratio, since the class grows only as ln(1 + scan*p)/p while the cost
    is linear in it.
    """
    free = torch.cuda.mem_get_info()[0]
    budget = frac * free - pchunk * block * 2
    per_line = 72 + 4 + (72 * 2 * 2) / float(ratio)          # pool + index + window and copy
    pool = int(min(cap, max(4_000_000, budget / per_line)))
    return pool, max(1_000_000, pool // ratio)


# ---------------------------------------------------------------- pool generation

def gen_seeds(nseed, seed=0, orbit=2_100_000, verbose=True):
    """Seeds spread over MANY automorphism orbits.

    A breadth-first closure of the 72 basis vectors under the generators produces ONE orbit,
    and one orbit is atypically self-conflicting: measured conflict density 0.010241 against
    0.005469 for a mixed set, which halves every class built from it.  The production rule
        u, v minimal with <u,v> = 4  =>  |u-v|^2 = 8 + 8 - 8 = 8, so u-v is minimal
    leaves the orbit, and iterating it reaches many orbits quickly.

    THE ORBIT CLOSURE MUST BE CUT SHORT.  A random word in the generators maps each seed inside
    its OWN orbit, so the orbit diversity of the whole pool is decided here and nowhere else: if
    the closure is allowed to supply all `nseed` seeds they are all in one orbit and the finished
    pool measures p = 0.0104 instead of 0.0055, which costs almost half of every class.  So the
    closure is stopped at `orbit` lines and the difference rule supplies the rest.
    """
    dev = device()
    G = load_gram(); GE = load_gens()
    assert np.abs(GE).max() <= 127, "generator entries do not fit int8"
    Gt = torch.tensor(G.T.astype(np.int8), device=dev)          # x -> x G, entries <= MAXG
    GEt = [torch.tensor(A.astype(np.int8), device=dev) for A in GE]
    cur = torch.tensor(canon(np.eye(72, dtype=np.int64)).astype(np.int8), device=dev)
    nm = np.einsum('ij,jk,ik->i', np.eye(72, dtype=np.int64), G, np.eye(72, dtype=np.int64))
    assert (nm == 8).all(), "basis vectors are not minimal"
    rng = np.random.default_rng(seed)
    r1 = rng.integers(1, 2 ** 40, size=72); r2 = rng.integers(1, 2 ** 40, size=72)
    Xout = torch.empty((nseed, 72), dtype=torch.int8, device=dev)
    Hs = torch.empty(0, dtype=torch.int64, device=dev)
    n, frontier = 0, cur
    n, Hs = dedup_into(cur, _hash_keys(cur, r1, r2), Xout, Hs, n)
    # phase 1: automorphism closure, phase 2: the difference rule on the frontier
    stop = min(orbit, nseed)
    while n < stop:
        grew = 0
        for A in GEt:
            if n >= stop:
                break
            V = canon_t(mm(frontier, A).to(torch.int8))
            m, Hs = dedup_into(V, _hash_keys(V, r1, r2), Xout, Hs, n)
            grew += m - n; n = m
        if grew == 0:
            break
        frontier = Xout[n - grew:n].clone()
        if verbose:
            print("   seeds (orbit) %d" % n, flush=True)
    while n < nseed:
        k = min(8192, n)
        idx = torch.tensor(np.sort(rng.choice(n, k, replace=False)), device=dev)
        S = Xout[idx]
        M = mm(S, mm(S, Gt).to(torch.int8).T.contiguous())
        ii, jj = torch.nonzero(torch.triu(M == 4, 1), as_tuple=True)
        if len(ii) == 0:
            break
        D = canon_t((S[ii].to(torch.int16) - S[jj].to(torch.int16)).to(torch.int8))
        m, Hs = dedup_into(D, _hash_keys(D, r1, r2), Xout, Hs, n)
        if m == n:
            break
        n = m
        if verbose:
            print("   seeds (difference rule) %d" % n, flush=True)
    return Xout[:n]


def diff_sweep(P, nu, rng, chunk=2_000_000, block=64, maxrows=24_000_000):
    """New minimal lines by the production rule, using WHOLE-POOL neighbourhoods.

        u, v minimal, <u,v> = 4  =>  |u-v|^2 = 8 + 8 - 8 = 8, so u-v is minimal
        u, v minimal, <u,v> = -4 =>  u+v is minimal, by the same count

    The point is which pairs are looked at.  Sampling k lines and taking the ip = +-4 pairs
    inside the sample finds only C(k,2) * P(|ip| = 4) of them, which for k = 8192 is about
    50 000 and dries up fast.  Sweeping ONE u against the whole pool finds n * P(|ip| = 4)
    neighbours at once -- about 30 000 for every single u at n = 20M -- and a block of 64 u's
    costs a single pass over the pool.  This is also the ONLY step that leaves an automorphism
    orbit, so it is what decides how much of the shell is reachable at all: a random word maps
    every seed inside its own orbit, and one orbit holds only about 2.06M of the 3.1e9 lines.
    """
    dev = P.device
    Gt = torch.tensor(load_gram().T.astype(np.int8), device=dev)
    if len(P) > maxrows:                    # sweeping the WHOLE pool makes the generator
        off = int(rng.integers(0, len(P)))  # superlinear in the target: the yield per u is
        sl = (off + np.arange(maxrows)) % len(P)          # proportional to the rows scanned,
        P = P[torch.tensor(sl, device=dev)]               # so cap the rows and raise nu
    out = []
    got = 0
    for a in range(0, nu, block):
        m = min(block, nu - a)
        ui = torch.tensor(np.sort(rng.choice(len(P), m, replace=False)), device=dev)
        U = P[ui]
        Y = mm(U, Gt).to(torch.int8).T.contiguous()               # (72, m), entries <= MAXG
        for c in range(0, len(P), chunk):
            ip = mm(P[c:c + chunk], Y)
            ii, jj = torch.nonzero(ip.abs() == 4, as_tuple=True)
            if not len(ii):
                continue
            sg = torch.sign(ip[ii, jj]).to(torch.int16)
            D = U[jj].to(torch.int16) * sg[:, None] - P[c + ii].to(torch.int16)
            out.append(canon_t((-D).to(torch.int8)))
            got += len(ii)
            del ip, ii, jj, sg, D
    return torch.cat(out) if out else P[:0]


def spread_pool(target, seed=0, orbit=2_100_000, walk=30, nu=512, sample=2_000_000,
                batch=8_000_000, maxrows=24_000_000, nseed=None, verbose=True):
    """A pool of `target` distinct minimal LINES spread over many automorphism orbits.

    TWO MOVES, AND EACH DOES SOMETHING THE OTHER CANNOT.

      walk   push a random sample of the pool by an independent random word in the generators.
             Uniform measure is stationary for a random walk on a finite group, so a long word
             lands close to uniformly in each line's own orbit; one word per batch (rather than
             per row) makes every step a single dense GEMM.  This FILLS orbits but can never
             leave one, and |Aut| bounds an orbit at about 2.06M of the 3.1e9 lines -- so the
             walk alone saturates, which is exactly what a first version did at ~25M.
      diff   the production rule <u,v> = +-4  =>  u -+ v minimal, swept over the whole pool.
             This is the only move that OPENS a new orbit, and each sweep pays for itself: at
             n = 20M every u has about 2900 neighbours in the pool, so 512 u's cost one pass
             and return 1.5M lines.

    They are alternated, and the walk sample is drawn from the CURRENT pool rather than from a
    fixed seed set, so orbits opened by a sweep are filled by the next walk.  The sweep reads
    at most `maxrows` of the pool, since its yield is proportional to the rows it scans and
    letting it grow with the pool makes the whole generator superlinear in the target.  The conflict
    density that results, p = 0.00546, is the mixed-shell figure; one orbit alone gives 0.0102
    and halves every class built from it.
    """
    dev = device()
    G = load_gram()
    GEt = [torch.tensor(A.astype(np.int8), device=dev) for A in load_gens()]
    rng = np.random.default_rng(seed)
    r1 = rng.integers(1, 2 ** 40, size=72); r2 = rng.integers(1, 2 ** 40, size=72)
    Xout = torch.empty((target, 72), dtype=torch.int8, device=dev)
    Hs = torch.empty(0, dtype=torch.int64, device=dev)
    S = gen_seeds(min(orbit, target), seed=seed, orbit=orbit, verbose=False)
    n, Hs = dedup_into(S, _hash_keys(S, r1, r2), Xout, Hs, 0)
    del S
    if verbose:
        print("   orbit closure: %d lines" % n, flush=True)
    stall = 0
    while n < target and stall < 3:
        before = n
        t0 = time.time()
        D = diff_sweep(Xout[:n], nu, rng, maxrows=maxrows)
        if len(D):
            n, Hs = dedup_into(D, _hash_keys(D, r1, r2), Xout, Hs, n)
        del D
        mid = n
        if verbose:
            print("   diff  +%-9d -> %d / %d  (%.1fs)"
                  % (n - before, n, target, time.time() - t0), flush=True)
        made = 0
        t0 = time.time()
        while made < batch and n < target:
            k = min(sample, n)
            idx = torch.tensor(np.sort(rng.choice(n, k, replace=False)), device=dev)
            W = Xout[idx]
            for g in rng.integers(0, len(GEt), size=walk):
                W = mm(W, GEt[g]).to(torch.int8)
            assert int(W.abs().max()) <= MAXC, "walk left the shell"
            b4 = n
            n, Hs = dedup_into(canon_t(W), _hash_keys(canon_t(W), r1, r2), Xout, Hs, n)
            made += k
            del W, idx
            if n == b4:
                break
        if verbose:
            print("   walk  +%-9d -> %d / %d  (%.1fs)"
                  % (n - mid, n, target, time.time() - t0), flush=True)
        stall = stall + 1 if n == before else 0
    X = Xout[:n]
    Gi = torch.tensor(G, dtype=torch.float64, device=dev)
    for a in range(0, n, 2_000_000):                    # exact re-check of every norm
        B = X[a:a + 2_000_000].to(torch.float64)
        assert bool((((B @ Gi) * B).sum(1) == 8).all()), "pool row of norm != 8"
    return X


# The inner-product distribution of the Gamma_72 shell is known EXACTLY.  Its minimal
# vectors form a spherical 11-design (Venkov, for an extremal even unimodular lattice), and
# the attainable products are 0, +-1, ..., +-4, +-8, so five moments plus the total pin the
# distribution down, with the tenth moment left over as a check.  Per minimal vector:
#     <u,v> = 0 : 2 603 658 750     +-1 : 1 512 243 200 each     +-2 : 280 928 256 each
#            +-3 :    13 959 168 each      +-4 :     127 800 each      +-8 : 1
# Identifying v with -v, the number of LINES at each |<u,v>| and the fractions are:
SHELL_LINES = {0: 1301829375, 1: 1512243200, 2: 280928256, 3: 13959168, 4: 127800}
SHELL_TOTAL = sum(SHELL_LINES.values())                       # 3 109 087 799
SHELL_FRAC = {c: n / SHELL_TOTAL for c, n in SHELL_LINES.items()}
SHELL_P = (SHELL_LINES[3] + SHELL_LINES[4]) / SHELL_TOTAL     # 0.0045309 exactly


def ip_histogram(X, nprobe=64, seed=0, chunk=4_000_000):
    """measured |<u,v>| distribution of the pool against the exact shell distribution.

    This is the pool's quality certificate and the one number that drives everything: greedy
    returns about ln(1 + n p)/p lines, so a pool 20% denser than the shell loses 20% of every
    class.  A uniform random sample of the shell would reproduce SHELL_FRAC; a single
    automorphism orbit gives p = 0.0102, more than twice the shell's 0.0045309.
    """
    dev = X.device
    Gt = torch.tensor(load_gram().T.astype(np.int8), device=dev)
    rng = np.random.default_rng(seed)
    idx = torch.tensor(np.sort(rng.choice(len(X), nprobe, replace=False)), device=dev)
    Y = mm(X[idx], Gt).to(torch.int8).T.contiguous()
    cnt = torch.zeros(9, dtype=torch.int64, device=dev)       # counted ON the device: the
    for a in range(0, len(X), chunk):                         # products are 23 GB at n = 90M
        v = mm(X[a:a + chunk], Y).abs().clamp_(max=8).reshape(-1).to(torch.int64)
        cnt += torch.bincount(v, minlength=9)[:9]
    cnt = cnt.cpu().numpy()
    cnt[8] -= nprobe                                          # each probe against itself
    tot = int(cnt.sum())
    return {c: int(cnt[c]) for c in range(9)}, tot


def conflict_density(X, nprobe=64, seed=0):
    """density of pairs with |<u,v>| >= 3 -- the p that fixes the class size"""
    h, tot = ip_histogram(X, nprobe=nprobe, seed=seed)
    return sum(h[c] for c in range(3, 9)) / float(tot)


# ---------------------------------------------------------------- disjoint classes

def _greedy_block(M, thr):
    """maximal independent subset of a small block, from its own Gram; plain sequential greedy.
    Blocking does not change the answer: the block is drawn from the alive set, so it is
    already compatible with everything accepted before, and it is made internally compatible
    here.  The result is exactly what one-at-a-time greedy in the same order would give."""
    b = len(M)
    alive = np.ones(b, bool)
    keep = []
    for i in range(b):
        if not alive[i]:
            continue
        keep.append(i)
        alive &= (np.abs(M[i]) <= thr)
        alive[i] = False
    return np.array(keep, dtype=np.int64)


def pool_digest(X, stride=79, chunk=2_000_000):
    """A fingerprint of the pool: SHA-256 of every `stride`-th row, in pool order.

    WHY THIS EXISTS.  Everything downstream is indexed BY POOL POSITION -- the alive bitmask a
    --resume restores, the per-class digests, the SHA-256 of the whole family.  None of that
    means anything unless the pool a later run regenerates from the seed is the pool the
    earlier one used.  The pipeline is deterministic within one software stack, but it is NOT
    guaranteed across them: the pool is assembled by torch primitives whose tie-breaking is
    unspecified, and a run resumed against a DIFFERENT pool would silently reuse lines that the
    bitmask says are free -- the invariant |alive| + sum(sizes) = |pool| counts rows and would
    not notice.  So the digest is recorded and a resume refuses to proceed without it.

    A strided sample rather than the whole 5.8 GB: one pass over 1/79 of the rows costs a
    second instead of a minute, and any reordering or substitution of pool rows changes it with
    overwhelming probability.  `len(X)` is checked separately and exactly."""
    h = hashlib.sha256()
    h.update(("gamma72-pool-v1 n=%d stride=%d " % (len(X), stride)).encode())
    for a in range(0, len(X), chunk * stride):
        blk = X[a:a + chunk * stride:stride]
        h.update(blk.cpu().numpy().tobytes())
    return h.hexdigest()


def build_classes(X, nclasses, seed=0, thr=2, scan=20_000_000, block=256, pchunk=2097152,
                  runlen=65536, scanbig=None, nbig=0, usefused=False, permute=False,
                  cachefree=50, state=None, resume=False,
                  verbose=True, report=50, checkpoint=None):
    """Pairwise-disjoint classes of Gamma_72 minimal lines, greedy, on the GPU.

    Per class a window of `scan` still-unused lines is gathered into a contiguous float16
    buffer and greedily reduced.  Two parameters and why they are what they are:

      scan   greedy returns about ln(1 + scan*p)/p lines, so the class grows only
             LOGARITHMICALLY in the window while the cost grows linearly.
      scanbig/nbig
             the first `nbig` classes use a window of `scanbig` instead.  Dimension 72+k reads
             only the T+P largest classes -- 102 of them for dimension 81, 31 072 for
             dimension 95 -- so making the first couple of thousand as large as possible is
             what the low dimensions live on, and it costs almost nothing: 2500 classes out of
             31 072 barely touch the pool, and the remaining classes are unaffected.
      runlen the window is scan/runlen CONTIGUOUS RUNS of the alive list taken at random
             offsets, not scan individual random rows.  The alive list is in pool order and
             the pool is written in hash order inside each generation batch, so a run is
             already a random subset of its batch and a few hundred runs mix across batches;
             but the gather now reads a few hundred contiguous stretches of the 6.5 GB pool
             instead of scan scattered rows, which on a T4 is the difference between 1.59 and
             0.5 seconds a class.  runlen = 1 recovers the fully scattered version.
      block  the round costs (product) 4*b*alive bytes and (operand + compaction) 288*alive
             bytes, and alive sums to scan*(1+pb)/(pb) over the rounds, so the total is
             4*scan*(1+pb)/p + 288*scan*(1+pb)/(pb): the first term wants b small, the second
             wants b large, and the sum is flattest near b = 128.
      pchunk rows per product tile.  The (alive x block) product is 25 GB of the 32 GB a class
             moves, and it is written and then immediately reduced away, so tiling it to sit
             inside the T4's 4 MB L2 (pchunk * block * 2 bytes) keeps it off HBM entirely.
             The cost is one extra kernel launch a tile, so there is an optimum.

    The operand is float16 and so is the product.  That is exact for the only test made: the
    inputs (|x| <= 24, |xG| <= 8) are integers well inside float16's exact range 2048, the
    accumulation is float32 (13 824 << 2^24), and a stored value of modulus >= 2048 comes back
    within 8 of the truth, so it is never mistaken for one of modulus <= 2.  assert_fp16_exact
    re-checks that against int64 before any class is built.

    An int8 operand halves the gather and compaction traffic but torch._int_mm returns int32,
    which DOUBLES the product traffic -- and the product dominates.  Measured on the T4 at
    scan 20M: float16 0.76 s a class, int8 1.57 s.  A fused triton kernel would remove the
    product from memory altogether; on sm_75 its int8 tl.dot does not compile at all and the
    float16 one, though exact, timed slower than cuBLAS, so it is off by default and
    fused.select() will only ever return it if it measures faster (see fused.py).
    """
    dev = X.device
    G = load_gram()
    assert_fp16_exact(dev)
    fused = _fused.select(dev, thr=thr, verbose=verbose) if usefused else None
    CH = pchunk
    Gh = torch.tensor(G.T, dtype=torch.float16, device=dev)       # x -> x @ Gh = x G
    Gt8 = torch.tensor(G.T.astype(np.int8), device=dev)
    rng = np.random.default_rng(seed)
    # The alive list is kept in POOL ORDER, and deleting preserves order.  It must not be
    # permuted: a run of it would then be `runlen` random rows of a 6 GB pool, which is the
    # scattered gather the run window exists to avoid.  Randomness comes from the run offsets
    # and from the pool's own order, which is the hash order each dedup batch is written in.
    live = (torch.tensor(rng.permutation(len(X)).astype(np.int32), device=dev)
            if permute else torch.arange(len(X), dtype=torch.int32, device=dev))
    sizes, digests, worst, c0 = [], [], 0, 0
    pdig = pool_digest(X)
    if resume and state and os.path.exists(state):
        # np.load on an npz is lazy and keeps the file open, which on Windows blocks the
        # os.replace that rotates the state later in the run.  Read everything, then close.
        with np.load(state, allow_pickle=True) as z:
            keep = np.unpackbits(z['alive'])[:len(X)].astype(bool)
            sizes = [int(v) for v in z['sizes']]
            digests = [bytes(d) for d in z['digests']]
            worst = int(z['worst'])
            rngstate = str(z['rng'])
            saved = str(z['pool_digest']) if 'pool_digest' in z.files else None
        # THE RESUME GUARD.  The bitmask indexes pool POSITIONS, so resuming against a pool
        # that is not bit-for-bit the one it was written against would reuse lines already
        # taken, and no invariant in this file would catch it.  A state written before this
        # check existed carries no digest, and is refused rather than trusted.
        if saved is None:
            raise SystemExit(
                "%s carries no pool digest, so it cannot be shown to belong to this pool.  "
                "It was written by a version of this file that did not record one; rebuild "
                "from scratch rather than resuming." % state)
        if saved != pdig:
            raise SystemExit(
                "POOL MISMATCH: the state was written against pool %s, this run regenerated "
                "%s.  The alive bitmask indexes pool positions, so resuming would silently "
                "reuse lines.  Rebuild from scratch." % (saved[:16], pdig[:16]))
        live = torch.tensor(np.nonzero(keep)[0].astype(np.int32), device=dev)
        del keep
        c0 = len(sizes)
        rng.bit_generator.state = json.loads(rngstate)
        if verbose:
            print("   resumed at class %d, %d lines still unused, worst |ip| so far %d"
                  % (c0, len(live), worst), flush=True)
    out = []
    t_last = time.time()
    for c in range(c0, nclasses):
        nl = len(live)
        if nl < 20000:
            break
        s = min(scanbig if (scanbig and c < nbig) else scan, nl)
        rl = min(runlen, s)
        nb = max(1, s // rl)
        offs = torch.tensor(rng.integers(0, nl, size=nb), dtype=torch.int64, device=dev)
        w = (offs[:, None] + torch.arange(rl, dtype=torch.int64, device=dev)[None, :]) % nl
        w = torch.unique(w.reshape(-1))              # runs may overlap; the window is a set
        Xf = X[live[w].to(torch.int64)].to(torch.float16)
        cand = w.to(torch.int32)
        acc = []
        while len(cand):
            bb = min(block, len(cand))
            Sb = Xf[:bb]
            Yb = (Sb @ Gh).T.contiguous()                        # (72, bb), |entries| <= MAXG
            M = (Sb @ Yb).to(torch.int32).cpu().numpy()
            keep = _greedy_block(M, thr)
            kt = torch.tensor(keep, dtype=torch.long, device=dev)
            acc.append(cand[kt])
            Yt = Yb[:, kt].contiguous()
            ok = torch.empty(len(Xf), dtype=torch.bool, device=dev)
            for a0 in range(0, len(Xf), CH):
                ok[a0:a0 + CH] = (Xf[a0:a0 + CH] @ Yt).abs().amax(dim=1) <= thr
            ok[:bb] = False
            Xf = Xf[ok]
            cand = cand[ok]
        sel = torch.cat(acc)
        C = torch.sort(live[sel.to(torch.int64)].to(torch.int64)).values
        keepmask = torch.ones(nl, dtype=torch.bool, device=dev)
        keepmask[sel.to(torch.int64)] = False
        live = live[keepmask]
        Cn = C.cpu().numpy()
        # DISJOINTNESS INVARIANT.  Every class is drawn from `live` and its members are then
        # struck out of it, so |live| + sum(sizes) = |pool| holds after each class if and only
        # if no line was ever taken twice.  Checking it here makes disjointness a running
        # invariant rather than something to re-derive at the end -- which matters because a
        # resumed run no longer holds the earlier classes' lines.
        assert len(live) + sum(sizes) + len(Cn) == len(X), "a line was used twice"
        # Verify NOW, in exact integer arithmetic, and keep only the size and a digest: a
        # 31 000-class family is 250 MB of lines, and none of it is needed afterwards.
        worst = max(worst, verify_one(X, Cn, Gt8))
        sizes.append(int(len(Cn)))
        digests.append(hashlib.sha256(np.sort(Cn).astype(np.int64).tobytes()).digest())
        out.append(Cn)
        del Xf, cand, acc, sel, keepmask
        # The working window is resized once a round -- 20-odd rounds a class, tens of
        # thousands of classes -- so the caching allocator ends up holding blocks of hundreds
        # of distinct sizes.  A 60-class run stays at 0.29 s a class; by class 500 the same
        # parameters had drifted to 6.9 s.  Handing the blocks back periodically costs one
        # synchronisation and a few cudaMallocs.
        if cachefree and (c + 1) % cachefree == 0 and dev.type == 'cuda':
            torch.cuda.empty_cache()
        if (c + 1) % report == 0:
            if verbose:
                now = time.time()
                print("   %6d classes, sizes %d..%d, %d lines, %d unused, %.3f s/class"
                      % (c + 1, min(sizes), max(sizes),
                         sum(sizes), len(live), (now - t_last) / report), flush=True)
                t_last = now
            if checkpoint is not None:
                checkpoint(sizes)        # a long run must survive losing its session
            if state:
                # Written to a temporary file and renamed, because the state is copied off the
                # machine while the run continues: a reader that catches a half-written npz
                # gets BadZipFile and the run has no recoverable checkpoint at all.  os.replace
                # is atomic within a filesystem, so a reader sees either the old file or the
                # new one.  The previous state is kept as .bak for the same reason.
                keep = np.zeros(len(X), dtype=bool)
                keep[live.cpu().numpy().astype(np.int64)] = True
                tmp = state + '.tmp'
                np.savez(tmp, alive=np.packbits(keep), sizes=np.array(sizes, dtype=np.int32),
                         digests=np.array(digests), worst=worst, pool_digest=pdig,
                         rng=json.dumps(rng.bit_generator.state, default=str))
                if not tmp.endswith('.npz'):
                    tmp += '.npz'                      # np.savez appends the suffix
                try:
                    if os.path.exists(state):
                        os.replace(state, state + '.bak')
                except OSError:
                    pass                               # a stale reader; the rename below still runs
                os.replace(tmp, state)
    return out, sizes, digests, worst, pdig


def verify_one(X, C, Gt, chunk=4096):
    """exact check of a single class: norms 8, |off-diagonal| <= 2.  Returns the worst
    |<u,v>| found, so a long run can verify as it goes and never keep the lines."""
    V = X[torch.tensor(C, dtype=torch.long, device=X.device)]
    VG = mm(V, Gt)
    assert int(VG.abs().max()) <= MAXG, "x G outside its proved bound"
    nm = (VG.to(torch.int16) * V.to(torch.int16)).sum(1)
    assert bool((nm == 8).all()), "class contains a vector of norm != 8"
    VGb = VG.to(torch.int8).contiguous()
    Vt = V.T.contiguous()
    worst = 0
    for a in range(0, len(V), chunk):
        ip = mm(VGb[a:a + chunk], Vt)
        r = torch.arange(len(ip), device=X.device)
        ip[r, r + a] = 0
        worst = max(worst, int(ip.abs().max()))
    return worst


def verify_classes(X, classes, chunk=4096, float64_sample=8):
    """Exact re-check of a family: every vector of norm 8, every class a 60-degree code
    (|<u,v>| <= 2 off the diagonal), and all lines distinct across all classes.

    The arithmetic is int8 x int8 -> int32 on the tensor cores, which is exact by definition
    (|<x,y>| <= 72 * 24 * 8 = 13 824 fits int32 with room to spare) and about thirty times
    faster than float64 on a T4, which is what makes checking all 31 072 classes rather than
    a sample affordable.  A few classes are re-done in float64 as a cross-check that the
    integer path is wired up correctly.
    """
    dev = X.device
    Gt = torch.tensor(load_gram().T.astype(np.int8), device=dev)
    G64 = torch.tensor(load_gram(), dtype=torch.float64, device=dev)
    worst = 0
    for ci, C in enumerate(classes):
        V = X[torch.tensor(C, dtype=torch.long, device=dev)]
        VG = mm(V, Gt)
        assert int(VG.abs().max()) <= MAXG, "x G outside its proved bound"
        nm = (VG.to(torch.int16) * V.to(torch.int16)).sum(1)
        assert bool((nm == 8).all()), "class contains a vector of norm != 8"
        VGb = VG.to(torch.int8).contiguous()
        for a in range(0, len(V), chunk):
            ip = mm(VGb[a:a + chunk], V.T.contiguous())
            r = torch.arange(len(ip), device=dev)
            ip[r, r + a] = 0
            worst = max(worst, int(ip.abs().max()))
        if ci < float64_sample:                      # cross-check the integer path
            W = V.to(torch.float64)
            ip = (W @ G64) @ W.T
            ip.fill_diagonal_(0)
            assert int(ip.abs().max()) == max(worst, 0) or int(ip.abs().max()) <= 2
            assert bool((((W @ G64) * W).sum(1) == 8).all())
    assert worst <= 2, "class is not a 60-degree code: worst |ip| = %d" % worst
    allidx = np.concatenate(classes)
    assert len(np.unique(allidx)) == len(allidx), "classes are not pairwise disjoint"
    return worst
