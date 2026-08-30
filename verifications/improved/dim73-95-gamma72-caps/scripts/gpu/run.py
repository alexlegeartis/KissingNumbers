#!/usr/bin/env python3
"""Build pairwise-disjoint Gamma_72 classes on a GPU and emit a COMPACT CERTIFICATE.

    python run.py --pool 80000000 --classes 31050 --scan 20000000 --block 128 --seed 0 --tag prod

WHAT IS BEING MAXIMISED.  Dimension 95's cap construction gains
    4 * (sum of the T largest sizes among pairwise-disjoint classes)   (T = 31 050 for k = 23)
so the deliverable is the class SIZES, and the lines themselves are only needed to prove those
sizes are attainable.  A class is a set of Gamma_72 minimal lines with pairwise |<u,v>| <= 2.

WHY A SEED IS A CERTIFICATE.  Every random choice is made by numpy's PCG64 from --seed and only
then moved to the device; the GPU does nothing but exact integer GEMM (int8 tensor cores, or
float16 in with float32 accumulate, whichever reproduces an int64 reference -- both are checked
before any class is built).  So the pipeline is deterministic and the 5 GB of lines need not
be shipped: regenerate.py rebuilds them from the seed and checks the hash.

    HOW FAR THAT GOES.  Every random DECISION is numpy's, but the pool is ASSEMBLED by torch
    primitives -- a sort by 64-bit key inside the deduplicator -- and the pool's row order is
    what the greedy window is drawn from.  Reproduction is therefore exact within one software
    stack and is not guaranteed across them; measured, a certificate made on a T4 under one
    Colab image did not reproduce under a later one.  Everything downstream is indexed by pool
    POSITION, so the certificate and the resumable state now carry a POOL DIGEST and a resume
    refuses to run against a pool that does not match it (gamma72.pool_digest).  Without that
    guard a resume against a different pool would silently reuse lines, and the invariant
    |alive| + sum(sizes) = |pool| counts rows and would not notice.

Writes to --out (default $OUT_DIR, else ./out):
    cert_<tag>.json      parameters, seed, |ip| histogram, class sizes, SHA-256, timings
    sizes_<tag>.npy      class sizes as int32          (about 120 kB for 31 200 classes)
    state_<tag>.npz      resumable state, about 11 MB  (--resume picks it up)

RUNS ANYWHERE.  Plain numpy and torch; no Colab, no notebook, nothing device-specific -- the
exact-GEMM back end is chosen by probing the device and checking it against int64, and there
is a float32 fallback that works on any GPU and on the CPU.  On a machine with more memory
than a T4, `--auto` sizes the pool and the window to it:

    pip install -r requirements.txt
    python run.py --auto --classes 31200 --seed 7 --tag prod       # one command
    python run.py --auto --classes 31200 --seed 7 --tag prod --resume   # after any interruption
    python checkpool.py out/cert_prod.json                         # independent pool check
    python regenerate.py out/cert_prod.json                        # rebuild and match the hash
"""
import argparse, hashlib, json, os, time
# Must be set BEFORE torch initialises CUDA.  The class builder resizes its working buffer
# once a round (23 rounds a class, 31 072 classes), so the caching allocator accumulates
# blocks of hundreds of distinct sizes and eventually starts calling cudaFree/cudaMalloc,
# which synchronises: a 40-class run stays at 0.39 s a class while a long one drifts to 6 s.
os.environ.setdefault('PYTORCH_CUDA_ALLOC_CONF', 'expandable_segments:True')
import numpy as np
import gamma72 as g7

ap = argparse.ArgumentParser()
ap.add_argument('--pool', type=int, default=20_000_000)
ap.add_argument('--nu', type=int, default=512, help='u vectors per difference sweep')
ap.add_argument('--maxrows', type=int, default=24_000_000, help='pool rows a sweep reads')
ap.add_argument('--sample', type=int, default=2_000_000, help='pool rows pushed by one random word')
ap.add_argument('--walk', type=int, default=30)
ap.add_argument('--orbit', type=int, default=2_100_000, help='seeds taken from the plain orbit closure; the rest come from the difference rule, which is the only source of orbit diversity')
ap.add_argument('--classes', type=int, default=31050)
ap.add_argument('--scan', type=int, default=20_000_000)
ap.add_argument('--block', type=int, default=256)
ap.add_argument('--pchunk', type=int, default=2097152)
ap.add_argument('--runlen', type=int, default=65536)
ap.add_argument('--scanbig', type=int, default=0, help='window for the first --nbig classes')
ap.add_argument('--nbig', type=int, default=0)
ap.add_argument('--permute', action='store_true', help='permute the alive list (scattered window)')
ap.add_argument('--usefused', action='store_true', help='try the triton kernel (only used if it measures faster)')
ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--report', type=int, default=250)
ap.add_argument('--cachefree', type=int, default=50, help='release cached GPU blocks every N classes')
ap.add_argument('--state', default='', help='resumable state file (about 11 MB); survives a lost session')
ap.add_argument('--resume', action='store_true')
ap.add_argument('--tag', default='t4')
ap.add_argument('--out', default='', help='output directory (default: $OUT_DIR, else ./out)')
ap.add_argument('--auto', action='store_true', help='size --pool and --scan from the free GPU memory')
# (there is no --verify: every class is verified exactly as it is built, which costs 5 ms
#  and means the lines never have to be kept -- see build_classes and verify_one)
a = ap.parse_args()
OUT = a.out or os.environ.get('OUT_DIR', 'out')
os.makedirs(OUT, exist_ok=True)
import torch
dev = g7.device()
name = torch.cuda.get_device_name(0) if dev.type == 'cuda' else 'cpu'
free = torch.cuda.mem_get_info()[1] / 2 ** 30 if dev.type == 'cuda' else 0
print("device: %s (%s, %.1f GB)" % (dev, name, free), flush=True)
print("TF32 off: matmul=%s cudnn=%s" % (torch.backends.cuda.matmul.allow_tf32,
                                        torch.backends.cudnn.allow_tf32), flush=True)
back = g7.select_backend(dev)
if a.auto and dev.type == 'cuda':
    a.pool, a.scan = g7.autosize(a.block, a.pchunk)
    print("--auto: pool %d, scan %d" % (a.pool, a.scan), flush=True)
need = g7.residency(a.pool, a.scan, a.block, a.pchunk)
print("projected GPU residency: %.1f GB (pool int8 %.1f + window %.1f + tile %.1f)"
      % (need[0], need[1], need[2], need[3]), flush=True)
if dev.type == 'cuda':
    free = torch.cuda.mem_get_info()[0] / 2 ** 30
    if need[0] > 0.72 * free:
        print("WARNING: %.1f GB of a free %.1f GB.  Above about 70%% the caching allocator "
              "starts evicting and a class slows by an order of magnitude; lower --scan."
              % (need[0], free), flush=True)
rep = {'args': vars(a), 'device': name, 'backend': back}

t = time.time()
P = g7.spread_pool(a.pool, walk=a.walk, seed=a.seed, orbit=a.orbit, nu=a.nu,
                   sample=a.sample, maxrows=a.maxrows)
rep['t_pool'] = time.time() - t
rep['pool'] = int(len(P))
print("pool %d lines in %.0fs" % (len(P), rep['t_pool']), flush=True)
t = time.time()
hist, htot = g7.ip_histogram(P, seed=a.seed)
print("   histogram in %.0fs" % (time.time() - t), flush=True)
rep['ip_histogram'] = hist
rep['ip_pairs'] = htot
rep['p_conflict'] = sum(hist[c] for c in range(3, 9)) / float(htot)
print("  |ip|   measured      exact shell     ratio")
for c in range(5):
    print("   %d   %.7f     %.7f     %.4f"
          % (c, hist[c] / htot, g7.SHELL_FRAC[c], (hist[c] / htot) / g7.SHELL_FRAC[c]), flush=True)
print("conflict density p = %.6f against the shell's exact %.7f (%+.1f%%)"
      % (rep['p_conflict'], g7.SHELL_P, 100 * (rep['p_conflict'] / g7.SHELL_P - 1)), flush=True)
print("  -> greedy predicts ln(1+scan*p)/p = %.0f lines a class"
      % (np.log(1 + min(a.scan, len(P)) * rep['p_conflict']) / rep['p_conflict']), flush=True)

del_ = torch.cuda.empty_cache() if dev.type == 'cuda' else None
def _ckpt(sizes):
    np.save(os.path.join(OUT, 'sizes_%s.npy' % a.tag), np.array(sizes, dtype=np.int32))
t = time.time()
STATE = a.state or os.path.join(OUT, 'state_%s.npz' % a.tag)
if a.resume and not os.path.exists(STATE):
    # A resumed run gets its state back by being re-uploaded with the bundle, so look next to
    # the script as well as in the output directory.
    alt = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'state_%s.npz' % a.tag)
    if os.path.exists(alt):
        import shutil
        shutil.copy(alt, STATE)
        print("resume: seeded %s from the uploaded bundle" % STATE, flush=True)
C, sizes, digests, worst, pdig = g7.build_classes(
    P, a.classes, seed=a.seed, scan=a.scan, block=a.block, pchunk=a.pchunk, runlen=a.runlen,
    state=STATE, resume=a.resume,
    scanbig=a.scanbig or None, nbig=a.nbig, usefused=a.usefused, permute=a.permute,
    cachefree=a.cachefree, report=a.report, checkpoint=_ckpt)
rep['t_classes'] = time.time() - t
sz = np.array(sizes, dtype=np.int32)
print("classes %d in %.0fs (%.3f s/class), sizes %d..%d, %d lines"
      % (len(sz), rep['t_classes'], rep['t_classes'] / max(1, len(sz)), sz.min(), sz.max(),
         sz.sum()), flush=True)

# Every class was verified exactly as it was built (norms 8, |off-diagonal| <= 2), so all that
# is left is that the classes are pairwise DISJOINT -- which, the pool being duplicate-free
# (checkpool.py), is disjointness of pool indices.
allidx = np.concatenate(C) if C else np.zeros(0, dtype=np.int64)
# Disjointness holds by the invariant |live| + sum(sizes) = |pool| asserted after every class
# (see build_classes); for the classes built in THIS session it is re-checked directly too.
if len(C):
    assert len(np.unique(allidx)) == len(allidx), "classes are not pairwise disjoint"
disj = ("all %d classes pairwise disjoint (running invariant; %d built this session "
        "re-checked directly)" % (len(sz), len(C)))
rep['pool_digest'] = pdig            # see gamma72.pool_digest: without this the SHA-256
rep['worst_ip'] = int(worst)         # below, and the alive bitmask, name nothing
rep['n_verified'] = len(sz)
print("VERIFIED exactly while building: %d classes, worst |off-diagonal| = %d (<= 2); %s"
      % (len(sz), worst, disj), flush=True)

h = hashlib.sha256()
for d in digests:
    h.update(d)
srt = np.sort(sz)[::-1].astype(np.int64)
rep.update({'n_classes': int(len(sz)), 'sizes_min': int(sz.min()), 'sizes_max': int(sz.max()),
            'total_lines': int(sz.sum()), 'sha256': h.hexdigest(),
            'sec_per_class': rep['t_classes'] / max(1, len(sz)),
            'top31050': int(srt[:31050].sum()), 'sizes': [int(x) for x in sz]})
print("sum of the top %d sizes: %d  ->  dimension 95 gain 4x = %d"
      % (min(31050, len(srt)), rep['top31050'], 4 * rep['top31050']), flush=True)
np.save(os.path.join(OUT, 'sizes_%s.npy' % a.tag), sz)
json.dump(rep, open(os.path.join(OUT, 'cert_%s.json' % a.tag), 'w'), indent=1)
print(json.dumps({k: v for k, v in rep.items() if k != 'sizes'}, indent=1), flush=True)
