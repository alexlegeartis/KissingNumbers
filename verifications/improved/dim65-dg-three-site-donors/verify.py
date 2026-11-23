#!/usr/bin/env python3
"""Replay the full retained atlas and certify the three-profile replacements."""
import json, runpy, sys
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from hashlib import sha256
from tempfile import TemporaryDirectory
from subprocess import run
assert __debug__
here=Path(__file__).resolve().parent
old=here.parent/'dim65-68-dg-retained-atlas'
pin={
 'verify.py':'7f2bca666ba844610662f46ef9a1e868f9cafc80385ae59bc423bcfc27890f99',
 'verify_atlas.cpp':'1dac94dd0e8c13244d18c3f97de14786741910d000026334c200c93e607f4fe7',
 'atlas.inc':'95a6b6344cc7a4af5dece3be36e2bb89a311ade757bc9e932334aa9382b5ed85',
 'dg.py':'eefd262c41ae6d93c3407cf41222b0231d40e1e2cf14573b0691c0fece8f0eba',
 'assignment.txt':'12b23c1814c6b4879195c12619b899a77d4a2f614f4d462c79ec65be06bbade7',
 'claims.json':'fd18842c427c0e16ca3be745f8b44f0e2d93583d8ed53a4340ef184444123f26'
}
assert all(sha256((old/k).read_bytes()).hexdigest()==v for k,v in pin.items())
sys.path.insert(0,str(old))
with redirect_stdout(StringIO()):
 E=runpy.run_path(str(old/'verify.py'))
assert (65,40414,180224) in E['rows']
assert E['out']['lower_bounds']['65']==422858578
assert len(E['A'])==len(set(E['A']))==128 and len(E['C'])==len(set(E['C']))==1024
spec=here/'triple_assignment.txt'
lines=spec.read_text().splitlines()
assert len(lines)==23 and lines[0]=='21'
assert lines[-1].split()==(old/'assignment.txt').read_text().splitlines()[-1].split()
with TemporaryDirectory() as t:
 T=Path(t)
 for name in ('cubic80.txt','hyperplane_heads.txt'):
  (T/name).write_bytes(E['sources'][name])
 exe=T/'triple_atlas'
 run(['g++','-O3','-std=c++20',str(here/'atlas_cross.cpp'),'-o',str(exe)],check=True)
 run([str(exe),str(T/'cubic80.txt'),str(spec),str(T/'hyperplane_heads.txt')],check=True)

from fractions import Fraction as F
from itertools import product
P=[(-F(1129,200),F(5,4),F(5,4),F(0)),(F(7,5),-F(11,2),F(32,25),F(107,100)),(F(7,5),F(32,25),-F(11,2),-F(107,100))]
N=[61+sum(y*y for y in p) for p in P]
D=[sum(abs(y-1)for y in p[:3])for p in P]
def gate(dot,n,m):return dot<=0 or 4*dot*dot<=n*m
worst={}
def check(name,d,n,m):
 assert gate(d,n,m),(name,d,n,m)
 b=float(d)/float(n*m)**.5
 if b>worst.get(name,0):worst[name]=b
for i,(p,n) in enumerate(zip(P,N)):
 assert n<=96 and sum(y<0 for y in p[:3])==1 and min(abs(y)for y in p[:3])>=1
 check('dense',32+D[i],n,64)
 check('four-site',sum(abs(y)for y in p[:3])+max(1,abs(p[3])),n,4)
 check('tail axes',abs(p[3]),n,1)
 for s in range(8):
  h=[p[j]for j in range(3)if (s>>j)&1]
  check('sparse',16-len(h)+sum(abs(y)for y in h)-2*int(any(y<0 for y in h)),n,16)
 for j,q in enumerate(P):
  m=N[j]
  if i!=j:check('same owner',61+sum(x*y for x,y in zip(p,q)),n,m)
  check('same cell distinct owners',16+sum(abs(x*y-1)for x,y in zip(p[:3],q[:3]))+p[3]*q[3],n,m)
  check('different cells',32+D[i]+D[j]+p[3]*q[3],n,m)
 # Atlas vectors have 32 body signs of magnitude sqrt(2), and tail
 # magnitude 2*sqrt(2). At distance four all marked signs agree, so
 # the large negative multiplier cannot cost anything; otherwise the
 # verified next allowable signed distance is at least six.
 for dot in (24+sum(max(0,y-1) for y in p[:3])+2*abs(p[3]),
             20+D[i]+2*abs(p[3])):
  assert 8*dot*dot<=n*72
bank=len(E['A'])*len(E['C']);assert bank==131072
count=E['out']['lower_bounds']['65']-32*bank+2*21*bank
assert count==424169298
out={'status':'PASS','lower_bounds':{'65':count},'net_gain_vs_atlas':count-422858578}
assert out==json.loads((here/'claims.json').read_text())
print(json.dumps(out,sort_keys=True))
