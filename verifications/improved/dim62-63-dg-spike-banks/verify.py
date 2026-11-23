#!/usr/bin/env python3
"""Punctured Delsarte–Goethals spikes over exactly rephased RM80/63D."""
from hashlib import sha256
from pathlib import Path
from collections import defaultdict,Counter
from contextlib import redirect_stdout
from io import StringIO
import gzip,base64,json,struct,runpy
from dg import code
assert __debug__
here=Path(__file__).resolve().parent
parent=here.parent
p=parent/'dim62-63-guarded-rephase'
pins={'verify.py':'29957c51f761edff380d832dece779fd9b26a44581a85103cb873f61e92c81b9',
'witness62.bin.gz.b64':'321e121bbc6c3834f35d63e79399b3302a197898988cfa75edb0724001beb7a7',
'witness63.bin.gz.b64':'78abe2f27b448e7c973e4a1b4bd5e798151bb704c48f0e9a407835524a91081e'}
assert all(sha256((p/f).read_bytes()).hexdigest()==s for f,s in pins.items())
with redirect_stdout(StringIO()):
    D=runpy.run_path(str(p/'verify.py'))
prior=D['res']['63']
assert prior==230159620
M=(0,20); mask=sum(1<<i for i in M); chosen=1<<20
rows=D['dense']; cube=D['D']['phase'];quadratic=D['D']['quadratic']
# Confirm 0/20 ordering, and match the inherited indexed marked phases.
assert D['D']['marks63']==M and D['parity']==(1,0)
assert len(D['D']['phase'])==80
assert set(cube)<=set(rows) and set(quadratic)<=set(rows)
# Replay every old bank guard with phases from the tiny joint witness.
J=json.loads((here/'sign63.json').read_text())
assert J['bad']==[]
bits=int(J['bits'],16); N=0; idx=0
Q=sorted(s for s in D['sparse'] if (s&mask).bit_count()<=1)
seen=set(Q)
G=defaultdict(list)
b=gzip.decompress(base64.b64decode((p/'witness63.bin.gz.b64').read_text()))
while idx<len(b):
    s,k=struct.unpack_from('<QB',b,idx);idx+=9
    assert s not in seen and (s&mask).bit_count()<=1
    seen.add(s)
    pos=[1<<j for j in range(64)if (s&~mask)>>j&1]
    for _ in range(k):
        i,t=struct.unpack_from('<HH',b,idx);idx+=4
        g=sum(pos[j] for j in range(len(pos)) if (t>>j)&1)
        assert g.bit_count()==8 and g&(s&Q[i])==g
        G[s].append((g,(bits>>N&1)^1));G[Q[i]].append((g,bits>>N&1))
        N+=1
assert N==J['num']==13045 and bits>>N==0
for s in seen:
    assert s.bit_count()==16
    assert all((s&w).bit_count()%2==0 for w in rows)
    mode=s&mask
    a=(int(mode!=1<<20))
    b=(chosen&mode).bit_count()%2
    assert a!=b  # even full-parent restriction; deleted mark sum fixes body parity.
    H={0:0}
    for g,t in [(s,a)]+G[s]:
        if g in H:assert H[g]==t
        else:H.update({x^g:y^t for x,y in list(H.items())})
    assert len(H) in (2,4,8,16,32)
    assert all(g in (0,s) or g.bit_count()==8 for g in H)
# Rephased constraints use the same physical overlap guards as the
# exact predecessor. Complete each balanced 16-column dual to rank 5;
# all used dual rows have weight 8 / 16, so dummy shortening halves
# the count precisely and keeps distance four.
C,A=code()
assert len(C)==1024 and len(A)==128
counts=Counter(a&mask for a in A)
assert counts=={0:32,1:32,1<<20:32,mask:32}
# DG quadratic parts and affine translations are in the actual RM2
# order (boolean truth table index in ordinary binary).
from functools import reduce
from operator import xor
R=set(quadratic)
assert all(sum(1<<i for i in range(64) if i&t==t) in R for t in range(64) if t.bit_count()<=2)
body=[i for i in range(64) if i not in M]
assert len(body)==62
c0=cube[:len(body)]
assert len(set(c0))==len(body)  # distinctly assigned parent cubic cosets
assert all(counts[(c^q^chosen)&mask]==32 for c in c0 for q in C)
def gate(d,a,b):return d<=0 or 4*d*d<=a*b
m=62-1+4**2
assert m==77
# All donor, predecessor, and late-layer pairs, before multiplying by radii.
assert gate(62-2*24+15,m,m)  # donor same spike (DG)
assert gate(62-2*16+6,m,m)   # donor different spike (different cubic)
assert gate(62-2*(16-1)+3,m,64)  # donor vs original dense equality
assert gate(16+3-2,m,16) and gate(15+3-2,m,16) # ordinary and one-mark banks
assert gate(4+3,m,4) and gate(2+3,m,4)  # ordinary and doubled-head affine planes
assert gate(4,m,1)  # body poles; donor's new tail is zero
gain=len(body)*len(C)*32
total=prior+gain
assert gain==2031616 and total==232191236
from verify62 import run as run62
r=run62(here,parent,p,C,A,cube)
assert {'status':'PASS','lower_bounds':{'62':r,'63':total}}==json.loads((here/'claims.json').read_text())
print(json.dumps(dict(status='PASS',lower_bound=total,gain=gain,old_guards=N,body_marks=len(body),signs_per_mark=len(C)*32)))
