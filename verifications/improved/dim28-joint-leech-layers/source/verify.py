"""Independent exact verifier for the Leech lifts exported by this research.

Uses Python unbounded integers in Q(sqrt(2),sqrt(3)). The only mathematical
input is the standard kissing property of the reconstructed Golay Leech shell.
It checks all changed geometric interactions, without numerical tolerances.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from collections import Counter
import itertools
import re

GEN = [
 [0,1,2,3,4,11,19,21], [0,1,2,3,5,15,20,23],
 [0,1,2,3,6,9,10,12], [0,1,2,3,7,13,18,22],
 [0,1,2,3,8,14,16,17], [0,1,2,4,5,9,13,17],
 [0,1,2,4,6,16,22,23], [0,1,2,4,7,8,12,20],
 [0,1,2,5,6,8,11,18], [0,1,2,5,7,10,16,21],
 [0,1,3,4,5,8,10,22], [0,2,3,4,5,6,7,14]]

def golay():
    words = np.zeros((1,24), dtype=np.int8)
    for support in GEN:
        row = np.zeros(24,dtype=np.int8)
        row[support] = 1
        words = np.concatenate([words, words ^ row])
    assert Counter(words.sum(1)) == {0:1,8:759,12:2576,16:759,24:1}
    return words

def shell():
    g=golay()
    a=[]
    for i,j in itertools.combinations(range(24),2):
        for x,y in itertools.product([-4,4],repeat=2):
            v=np.zeros(24,dtype=np.int8); v[i]=x; v[j]=y; a.append(v)
    b=[]
    signs=np.array(list(itertools.product([-1,1],repeat=8)),dtype=np.int8)
    signs=signs[np.prod(signs,axis=1)==1]
    for word in g[g.sum(1)==8]:
        v=np.zeros((128,24),dtype=np.int8)
        v[:,np.flatnonzero(word)]=2*signs; b.append(v)
    c=[]
    for i in range(24):
        v=1-2*g.copy(); v[:,i]*=-3; c.append(v)
    v=np.concatenate([np.array(a),*b,*c])
    assert v.shape==(196560,24)
    assert np.all((v.astype(np.int64)**2).sum(1)==32)
    return v

def parse_scalar(s):
    out=[0,0,0,0]
    for term in s.replace('-','+-').split('+'):
        if not term: continue
        m=re.fullmatch(r'(-?\d*)\*?sqrt\((2|3|6)\)',term)
        if m:
            a=m[1]; a=-1 if a=='-' else (1 if a=='' else int(a))
            out[{2:1,3:2,6:3}[int(m[2])]]+=a
        else: out[0]+=int(term)
    return out


ZERO=(0,0,0,0)
def add(x,y): return tuple(int(a)+int(b) for a,b in zip(x,y))
def sub(x,y): return tuple(int(a)-int(b) for a,b in zip(x,y))
def scale(x,n): return tuple(int(a)*int(n) for a in x)
def mul(x,y):
    a,b,c,d=map(int,x); e,f,g,h=map(int,y)
    return (a*e+2*b*f+3*c*g+6*d*h,a*f+b*e+3*c*h+3*d*g,a*g+c*e+2*b*h+2*d*f,a*h+d*e+b*g+c*f)
def sign2(a,b):
    if a==0: return (b>0)-(b<0)
    if b==0: return (a>0)-(a<0)
    if (a>0)==(b>0): return (a>0)-(a<0)
    t=a*a-2*b*b
    return ((a>0)-(a<0))*((t>0)-(t<0))
def sign(x):
    a,b,c,d=map(int,x); p=sign2(a,b); q=sign2(c,d)
    if p==0: return q
    if q==0 or p==q: return p
    delta=sign2(a*a+2*b*b-3*c*c-6*d*d,2*a*b-6*c*d)
    return p*delta
def dot(x,y):
    out=ZERO
    for a,b in zip(x,y): out=add(out,mul(a,b))
    return out
def compatible(s,n,m): return sign(s)<=0 or sign(sub(scale(mul(s,s),4),mul(n,m)))<=0

def parse(path):
    tails={}; classes=[]; extras=[]; expected=None; dimension=None
    with path.open(encoding='utf8') as f:
        for line in f:
            line=line.strip()
            if not line or line.startswith('#'): continue
            parts=line.split()
            if parts[0]=='dimension': dimension=int(parts[1])
            elif parts[0]=='points': expected=int(parts[1])
            elif parts[0]=='tail': tails[int(parts[1])]=(int(parts[2]),np.array([parse_scalar(s) for s in parts[3].split(',')],np.int64))
            elif parts[0]=='class':
                ids=list(map(int,parts[1].split(','))); n=int(parts[2])
                heads=np.array([list(map(int,next(f).strip().split(','))) for _ in range(n)],np.int64)
                classes.append((ids,heads))
            elif parts[0]=='pole': extras.append(np.vstack([np.zeros((24,4),np.int64),np.array([parse_scalar(s) for s in parts[1].split(',')],np.int64)]))
            elif parts[0]=='extra': extras.append(np.array([parse_scalar(s) for s in parts[1].split(',')],np.int64))
            else: raise ValueError(parts[0])
    return dimension,expected,tails,classes,extras

def verify(path):
    start=time.time(); dim,expected,tails,classes,extras=parse(path)
    assert classes[0][0]==[0] and np.all(tails[0][1]==0) and tails[0][0]==1
    m=shell(); keys={r.tobytes() for r in m}; allkeys=[]
    for ids,h in classes:
        assert h.shape[1]==24 and np.all((h*h).sum(1)==32)
        allkeys.extend(r.astype(np.int8).tobytes() for r in h)
        if ids!=[0]:
            g=h@h.T; np.fill_diagonal(g,-32); assert g.max()<=8
    assert len(allkeys)==196560 and len(set(allkeys))==196560 and set(allkeys)==keys
    print('PASS: full Leech-shell partition and every head-set angle',flush=True)
    # All normalized cap tails have norm one after division by 4d.
    for k,(d,t) in tails.items():
        assert len(t)==dim-24 and d>0
        if k: assert dot(t,t)==(16*d*d,0,0,0)
    tail_owner={tid:ci for ci,(ids,h) in enumerate(classes) for tid in ids}
    assert len(tail_owner)==len(tails)
    used=[k for k in tails if k]
    for ii,i in enumerate(used):
        di,ti=tails[i]
        for j in used[ii+1:]:
            dj,tj=tails[j]; ip=dot(ti,tj)
            if tail_owner[i]==tail_owner[j]: assert sign(add(ip,(8*di*dj,0,0,0)))<=0
            else: assert sign(sub(ip,(8*di*dj,0,0,0)))<=0
    print('PASS: every cap-tail inequality, hence all body/cap and cap/cap pairs',flush=True)
    body=classes[0][1]; norms=[dot(p,p) for p in extras]
    checks=0
    for e,(p,norm) in enumerate(zip(extras,norms)):
        assert p.shape==(dim,4) and sign(norm)>0
        active=np.flatnonzero(np.any(p[:24]!=0,axis=0))
        assert len(active)<=1, 'head must be an integer vector times one positive radical'
        basis=int(active[0]) if len(active) else 0
        h=p[:24,basis]
        # These matrix products use int64; explicit small bounds prevent overflow.
        assert int(np.abs(h).max())*24*4<2**62
        bh=int((body@h).max()); s=[0]*4; s[basis]=bh
        assert compatible(s,norm,(32,0,0,0)), ('equator',e)
        for ids,c in classes[1:]:
            ch=int((c@h).max())
            for tid in ids:
                d,t=tails[tid]; s=[0]*4; s[basis]=ch*d; s=add(s,dot(p[24:],t))
                assert compatible(s,norm,(48*d*d,0,0,0)), ('cap',e,tid)
                checks+=len(c)
        for j in range(e):
            assert compatible(dot(p,extras[j]),norm,norms[j]), ('auxiliary',e,j)
        checks+=len(body)+e
    total=sum(len(ids)*len(c) for ids,c in classes)+len(extras)
    assert total==expected
    result=dict(status='EXACT VERIFICATION PASSED',dimension=dim,points=total,auxiliary=len(extras),head_sizes=[len(c) for _,c in classes[1:]],explicit_interactions_covered=checks,arithmetic='unbounded integers in Q(sqrt(2),sqrt(3)); no tolerance',mathematical_input='standard kissing property of the generated Golay Leech minimal shell',sha256=hashlib.sha256(path.read_bytes()).hexdigest(),elapsed_seconds=round(time.time()-start,3))
    print(json.dumps(result,indent=2),flush=True)
    path.with_suffix('.verification.json').write_text(json.dumps(result,indent=2))
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('certificate',type=Path)
    verify(ap.parse_args().certificate)
