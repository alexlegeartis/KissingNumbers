"""Independently reconstruct all explicit signed vectors and replay in scalar C++17."""
from pathlib import Path
import gzip,json,struct,subprocess,tempfile
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from prep import input66
from itertools import combinations

def sign_code():
 gray=['00','10','11','01'] # low bit first
 rows=set()
 for t in ('01112','03110','21310','21132'):
  for m in (1,-1):
   q=''.join(gray[(m*int(c))%4]for c in t)
   for i in range(0,10,2):
    b=q[i:]+q[:i]
    rows.add(tuple(int(c)for c in b[1:]))
 assert len(rows)==40 and min(sum(a!=b for a,b in zip(x,y))for x,y in combinations(rows,2))==3
 return rows

def main():
 assert __debug__
 W,T,prior=input66()
 base=prior+json.loads((HERE/'supports.json').read_text())
 partial=json.loads(gzip.decompress((HERE/'partial.json.gz').read_bytes()))
 phases=sign_code()
 M=set(T)
 with tempfile.TemporaryDirectory(prefix='66D_scalar_audit_') as raw:
  d=Path(raw); inp=d/'witness.bin';exe=d/'audit'
  with inp.open('wb') as out:
   wr=lambda a:out.write(struct.pack('<Q',int(a)))
   mask=lambda a:(wr(a&((1<<64)-1)),wr(a>>64))
   wr(len(W));wr(len(T));wr(len(base)+len(partial))
   for a in W+T:mask(a)
   for s in base:
    p=[j for j in range(66)if s>>j&1]
    hits=[c for c in combinations(p,4)if sum(1<<j for j in c) in M]
    assert len(hits)<=1
    if hits:
     t=hits[0];u=[j for j in p if j not in t];p=[u[0],t[0],t[1],u[1],u[2],t[2],t[3],u[3],u[4]]
    mask(s);wr(40)
    for w in phases:mask(sum(1<<j for j,b in zip(p,w)if b))
   for s,ph in partial:
    p=[j for j in range(66)if s>>j&1];mask(s);wr(len(ph))
    for w in ph:mask(sum(1<<j for i,j in enumerate(p)if w>>i&1))
  subprocess.run(['g++','-O3','-march=native','-std=c++17',str(HERE/'audit.cpp'),'-o',str(exe)],check=True)
  r=json.loads(subprocess.check_output([str(exe),str(inp)],text=True))
 assert r['total']>342297172
 return r

if __name__=='__main__':print(json.dumps(main(),sort_keys=True))

