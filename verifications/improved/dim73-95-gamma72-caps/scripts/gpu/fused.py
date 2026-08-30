"""Fused |<x,y>| <= thr test: one pass over the alive rows, no (alive x block) product.

The class builder's cost is dominated by a matrix it never wants.  For a block of `b` accepted
lines it forms the (alive x b) table of inner products, reduces it to one boolean a row, and
throws it away; at b = 128 that table is 25 GB of the 46 GB a class moves on a T4.  This
kernel keeps the tile in registers and writes only the boolean, so a round costs 72 bytes a row
(the int8 operand) plus one, instead of 72 + 4b.

The operand is READ as int8 (72 bytes a row) and converted to float16 in registers, because
Turing's int8 `tl.dot` does not compile in triton (TritonGPUAccelerateMatmul fails on sm_75)
while its float16 MMA does.  That is still exact: the bound proved in gamma72.py is
|<x,y>| <= 72 * 24 * 8 = 13 824, the inputs (|x| <= 24, |xG| <= 8) are integers well inside
float16's exact range 2048, and the accumulation is float32, exact to 2^24.  select() checks
the kernel against an int64 reference on random data of the full attainable range and falls
back to torch if triton is missing, will not compile, or disagrees.
"""
import numpy as np, torch

HAVE = False
try:
    import triton
    import triton.language as tl
    HAVE = True
except Exception:                                            # no triton on this machine
    pass

if HAVE:
    @triton.jit
    def _maxabs(X, Y, OUT, M, B, THR,
                K: tl.constexpr, BM: tl.constexpr, BB: tl.constexpr, KP: tl.constexpr):
        pid = tl.program_id(0)
        rm = pid * BM + tl.arange(0, BM)
        rk = tl.arange(0, KP)
        mk = rk < K
        mm_ = rm < M
        x = tl.load(X + rm[:, None] * K + rk[None, :],
                    mask=mm_[:, None] & mk[None, :], other=0).to(tl.float16)
        best = tl.zeros((BM,), dtype=tl.float32)
        for b0 in range(0, B, BB):
            rb = b0 + tl.arange(0, BB)
            y = tl.load(Y + rk[:, None] * B + rb[None, :],
                        mask=mk[:, None] & (rb[None, :] < B), other=0).to(tl.float16)
            acc = tl.dot(x, y, out_dtype=tl.float32)
            best = tl.maximum(best, tl.max(tl.abs(acc), axis=1))
        tl.store(OUT + rm, (best <= THR).to(tl.int8), mask=mm_)


CFG = (64, 64, 4)          # BM, BB, num_warps -- replaced by whatever select() times fastest


def ok_mask(X8, Y8, thr, cfg=None):
    """X8 (M,72) int8 rows, Y8 (72,B) int8 columns -> bool (M,), True iff max |ip| <= thr"""
    BM, BB, nw = cfg or CFG
    M, B = X8.shape[0], Y8.shape[1]
    out = torch.empty(M, dtype=torch.int8, device=X8.device)
    grid = ((M + BM - 1) // BM,)
    _maxabs[grid](X8, Y8, out, M, B, float(thr), K=72, BM=BM, BB=min(BB, max(16, B)),
                  KP=128, num_warps=nw)
    return out.bool()


CONFIGS = [(64, 64, 4), (128, 64, 4), (64, 32, 2), (128, 32, 4), (256, 64, 8), (64, 64, 8)]


def _time(fn, n=5):
    torch.cuda.synchronize()
    fn(); torch.cuda.synchronize()                            # warm up / compile
    t0 = torch.cuda.Event(True); t1 = torch.cuda.Event(True)
    t0.record()
    for _ in range(n):
        fn()
    t1.record(); torch.cuda.synchronize()
    return t0.elapsed_time(t1) / n


def select(dev, thr=2, verbose=True, bench_rows=4_000_000, bench_b=256, torch_mm=None):
    """Return ok_mask only if it is both EXACT and actually faster than the torch path.

    Being exact is not enough.  A first version compiled, verified, and then ran the class
    builder several times slower than plain int8 GEMM plus a reduction -- so the choice is
    made by measurement on a realistic shape, and the winning tile is picked the same way.
    """
    global CFG
    if not HAVE or dev.type != 'cuda':
        if verbose:
            print("   fused kernel: unavailable, using torch", flush=True)
        return None
    rng = np.random.default_rng(3)
    try:
        for M, B in ((1000, 128), (4097, 256), (77, 64), (33, 16)):
            A = rng.integers(-24, 25, size=(M, 72)).astype(np.int8)
            Bm = rng.integers(-8, 9, size=(72, B)).astype(np.int8)
            ref = np.abs(A.astype(np.int64) @ Bm.astype(np.int64)).max(axis=1) <= thr
            got = ok_mask(torch.tensor(A, device=dev),
                          torch.tensor(np.ascontiguousarray(Bm), device=dev), thr).cpu().numpy()
            assert (got == ref).all(), "fused kernel disagrees with the int64 reference"
    except Exception as e:
        if verbose:
            print("   fused kernel: rejected (%s: %s), using torch" % (type(e).__name__, e),
                  flush=True)
        return None
    X = torch.tensor(rng.integers(-24, 25, size=(bench_rows, 72)).astype(np.int8), device=dev)
    Y = torch.tensor(np.ascontiguousarray(rng.integers(-8, 9, size=(72, bench_b)).astype(np.int8)),
                     device=dev)
    best, bcfg = None, None
    for cfg in CONFIGS:
        try:
            t = _time(lambda: ok_mask(X, Y, thr, cfg))
        except Exception:
            continue
        if best is None or t < best:
            best, bcfg = t, cfg
    ref_t = _time(lambda: torch_mm(X, Y).abs().amax(dim=1) <= thr) if torch_mm else float('inf')
    del X, Y
    torch.cuda.empty_cache()
    if best is None or best >= ref_t:
        if verbose:
            print("   fused kernel: %s, using torch (%.1f ms against %.1f ms on %d x %d)"
                  % ("slower" if best else "no working tile", best or -1, ref_t, bench_rows,
                     bench_b), flush=True)
        return None
    CFG = bcfg
    if verbose:
        print("   fused kernel: active, tile %s, %.1f ms against torch's %.1f ms on %d x %d"
              % (bcfg, best, ref_t, bench_rows, bench_b), flush=True)
    return ok_mask
