#!/usr/bin/env python3
"""K(18) >= 8358: exact verification, self-contained.

    python verify18.py          # ~10 s; exits non-zero on any failure

Builds the configuration from data/ alone and checks every pair exactly.

Coordinates.  R^18 = R^16 (+) R^2.  Every vector has squared norm 8 (norm-8 units), so two
vectors are compatible iff their inner product is at most 4.  In the R^2 part every point has
one of the radii sqrt(8/9), sqrt2, sqrt8 and one of the angles 30k degrees, so every inner
product is (integer + integer*sqrt3)/36 after scaling by 36 (v-parts times 6, and the
radius-radius products 48, 72, 32, 144, 96, 288 times cos(30k)); the test a + b*sqrt3 <= 144
is decided by squaring.  No floating point is used in any decision.

The four families of vectors (counts 4320 + 3072 + 6 + 960 = 8358):
  equator   (v, 0):        v = (+-2,+-2,0^14) [480] and +-1 on each of the 30 octads of RM(1,4)
                            with an ODD number of minus signs [30 x 128]            (odd BW16)
  tiers     (a, p):        |p|^2 = 2 at angles 0,60,...,300; a = +-1 on a 6-set of the family
                            attached to that angle, odd number of minus signs [6 x 16 x 32]
  poles     (0, p):        |p|^2 = 8 at angles 0,60,...,300                                [6]
  tier B    ((2/3)(-1)^c, p): |p|^2 = 8/9 at angles 30,90,...,330, c a 16-bit word from
                            data/tierB.json                                          [6 x 160]
"""
import json, sys, os, itertools
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); DATA=os.path.join(HERE,'data')
def fail(msg): print('FAIL:',msg); sys.exit(1)
oct=json.load(open(os.path.join(DATA,'octads.json')))['octads']
fam=json.load(open(os.path.join(DATA,'families.json')))['families_at_angles_0_60_120']
tb=json.load(open(os.path.join(DATA,'tierB.json')))['slots_deg_to_words']
# ---- 1. the octads are RM(1,4): 30 words of weight 8 whose span with the all-ones word has dimension 5
def rank2(rows):
    M=np.array(rows,dtype=np.uint8)%2; r=0
    for c in range(M.shape[1]):
        p=[i for i in range(r,len(M)) if M[i,c]]
        if not p: continue
        M[[r,p[0]]]=M[[p[0],r]]
        for i in range(len(M)):
            if i!=r and M[i,c]: M[i]^=M[r]
        r+=1
    return r
def ind(s):
    v=np.zeros(16,dtype=np.uint8); v[list(s)]=1; return v
if len(oct)!=30 or any(len(o)!=8 for o in oct): fail('octads')
if rank2([ind(o) for o in oct])!=5: fail('octads do not span RM(1,4)')
if len(set(map(tuple,oct)))!=30: fail('octads repeated')
# ---- 2. build the vectors.  Scale v-parts by 6.  R^2 part = (r2, k): radius^2 in norm-8 units times 36, angle 30k.
groups=[]   # (int array n x 16, r2*36, k)
def signed(support, parity, scale):
    out=[]
    for bits in range(1<<len(support)):
        signs=[1-2*((bits>>t)&1) for t in range(len(support))]
        if sum(s<0 for s in signs)%2!=parity: continue
        v=np.zeros(16,dtype=np.int64)
        for t,i in enumerate(support): v[i]=scale*signs[t]
        out.append(v)
    return np.array(out)
eq=[]
for i,j in itertools.combinations(range(16),2):
    for si in (1,-1):
        for sj in (1,-1):
            v=np.zeros(16,dtype=np.int64); v[i]=12*si; v[j]=12*sj; eq.append(v)
eq=np.array(eq); eq=np.concatenate([eq]+[signed(o,1,6) for o in oct])
groups.append((eq,0,0))
for pos in range(6):
    F=fam[pos%3]
    if len(F)!=16 or any(len(s)!=6 for s in F): fail('family shape')
    groups.append((np.concatenate([signed(s,1,6) for s in F]),72,2*pos))
    groups.append((np.zeros((1,16),dtype=np.int64),288,2*pos))
for deg,words in tb.items():
    k=int(deg)//30
    V=np.array([[4*(1-2*((w>>(15-i))&1)) for i in range(16)] for w in words],dtype=np.int64)
    groups.append((V,32,k))
N=sum(len(g[0]) for g in groups)
# norms
for V,r2,k in groups:
    if not (((V*V).sum(1)+r2)==288).all(): fail('norm')
# distinct
keys=set()
for V,r2,k in groups:
    for v in V.tolist():
        t=(tuple(v),r2,k)
        if t in keys: fail('duplicate vector')
        keys.add(t)
# ---- 3. exact pair check.  cos(30m) = (alpha + beta*sqrt3)/2 with the table below.
COS={0:(2,0),1:(0,1),2:(1,0),3:(0,0),4:(-1,0),5:(0,-1),6:(-2,0),7:(0,-1),8:(-1,0),9:(0,0),10:(1,0),11:(0,1)}
bad=0; tight=0
for i,(A,ra,ka) in enumerate(groups):
    for j,(B,rb,kb) in enumerate(groups):
        if j<i: continue
        G=(A@B.T).astype(np.int64)
        if i==j: G=G[np.triu_indices(len(A),1)]
        R=int(round((ra*rb)**0.5))
        if R*R!=ra*rb: fail('radius product not a square')
        al,be=COS[(ka-kb)%12]
        # inner product = G + R*(al + be*sqrt3)/2 ;  need <= 144.  R is even in every case with be != 0.
        a2=2*G+R*al                    # twice the rational part
        b2=R*be                        # twice the sqrt3 coefficient (integer)
        # condition: a2 + b2*sqrt3 <= 288
        if b2>=0:
            v=(a2>288)|((a2<=288)&(3*b2*b2>(288-a2)**2))
        else:
            v=(a2>288)&(3*b2*b2<(a2-288)**2)
        bad+=int(np.sum(v)); tight+=int(np.sum((b2==0)&(a2==288)))
print(f'N = {N}   violating pairs = {bad}   pairs exactly at cos 1/2 (rational) = {tight}')
if bad: fail('kissing condition violated')
if N!=8358: fail('count')
# ---- 4. falsifiability: flipping one sign of one tier-B word must create violations
V,r2,k=groups[-1]; V2=V.copy(); V2[0,0]*=-1
G=(V2[:1]@eq.T).astype(np.int64)
if int(np.sum(G>144))==0: fail('self-test: perturbed word not caught')
print('OK: K(18) >= 8358, verified exactly (self-test passed)')
