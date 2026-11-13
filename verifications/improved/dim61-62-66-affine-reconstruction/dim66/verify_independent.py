from pathlib import Path
import json, hashlib, sys, itertools
import numpy as np
HERE=Path(__file__).resolve().parent
if not __debug__:raise RuntimeError('Run without -O; assertions are essential.')
ROOT=HERE.parents[3]; PKG=ROOT/'verifications/improved/dim61-67-compatible-sign-families'; SRC=HERE/'certificate.json'
c=json.loads(SRC.read_text());fam=c['family']; baseline=json.loads((PKG/'data/expanded66.json').read_text());assert fam==baseline
z=np.load(ROOT/'verifications/improved/dim61-68-95-96-code-layers/source/dim61_96/code64_supports.npz')
removed=0
z2=np.fromfile(ROOT/'verifications/improved/dim61-68-95-96-code-layers/source/dim61_96/cyclic73.weight16.bin',dtype='<u8').reshape(-1,2)
choice=json.loads((ROOT/'verifications/improved/dim61-68-95-96-code-layers/source/dim61_96/shortening_support_search.json').read_text())['dimensions']['66']
dual9=(21,37,38,49,57,63,67,70,72);coords=[i for i in range(1,73) if i not in dual9]+[0]+[i for i in dual9 if i not in choice['deleted_coordinates']]
parents=[]
for lo,hi in z2:
 w=int(lo)|(int(hi)<<64)
 if any(w>>i&1 for i in choice['deleted_coordinates']):continue
 parents.append(sum((w>>p&1)<<i for i,p in enumerate(coords)))
assert len(parents)==len(set(parents))==34028
extra=set(parents)-set(map(int,z['words']))
assert extra==set(json.loads((HERE/'context.json').read_text())['extra_full'])
G=[sum(int(x)<<i for i,x in enumerate(row)) for row in z['generator']]
cache={}
def build(r,full=False):
 s=r['support']; hs=r['checks']; t=r['shift']; assert isinstance(s,int) and 0<s<2**66 and s.bit_count()==16 and s&removed==0
 assert isinstance(t,int) and 0<=t<65536
 assert all(isinstance(h,int) and 0<=h<65536 for h in hs)
 ker=[w for w in range(65536) if all((w&h).bit_count()%2==0 for h in hs)]
 assert ker[0]==0 and min(w.bit_count() for w in ker if w)>=4
 assert 65535 in ker and all(w.bit_count()%2==0 for w in ker)
 assert len(ker)==(2048 if full else r['count'])
 if not full:assert len(ker)==2**r['dimension']
 coords=[i for i in range(66) if s>>i&1]
 words=[sum(((w^t)>>j&1)<<p for j,p in enumerate(coords)) for w in ker]
 assert len(set(words))==len(words)
 # Dense amplitudes alone bound every norm64 sparse/dense dot by32 in66D
 return {'s':s,'words':words,'n':len(words)}
fixed={r['support']:r for r in fam['old']+fam['new']};assert len(fixed)==len(fam['old'])+len(fam['new'])
mod={r['support']:r for r in c['modified_parent_classes']};assert len(mod)==len(c['modified_parent_classes']) and set(mod)<=set(parents) and not set(mod)&set(fixed)
added=[build(r) for r in c['new_classes']];assert len({a['s'] for a in added})==len(added) and not {a['s'] for a in added}&(set(parents)|set(fixed))
for s,r in mod.items():cache[s]=build(r,True)
for s in mod:
 assert all(s==t or (s&t).bit_count()<=8 for t in parents)
 assert all((s&r['support']).bit_count()<=8 for r in fam['new'])
# Independent direct sign-pair enumeration on all large intersections. Small ones
# satisfy dot <=4*8 regardless of sign.
count=critical=signpairs=0;maxdot=-64
oldlist=parents+[r['support'] for r in fam['new']]
for i,a in enumerate(added):
 for s in oldlist+[b['s'] for b in added[:i]]:
  count+=1; inter=a['s']&s; overlap=inter.bit_count()
  if overlap<=8:continue
  critical+=1
  if s not in cache:
   r=mod.get(s,fixed.get(s))
   if r is not None:cache[s]=build(r,True)
   else:
    cache[s]=next(b for b in added[:i] if b['s']==s)
  b=cache[s]
  # Reduce to distinct restricted sign patterns, preserving exact extrema.
  aa=set(w&inter for w in a['words']);bb=set(w&inter for w in b['words'])
  md=max(4*(overlap-2*(x^y).bit_count()) for x in aa for y in bb)
  assert md<=32,(a['s'],s,md)
  signpairs+=len(aa)*len(bb);maxdot=max(maxdot,md)
assert sum(a['n'] for a in added)==c['gain']==20288
# Each sparse point has16 coordinates +/-2, norm64; disjoint supports ensure
# no cross-class collisions. Dense points have60 nonzero body coordinates;
# quadruple/axis layers have4/1, so no layer collisions.
report=dict(status='PASS',gain=20288,totals={'66':342271364+20288},new_classes=len(added),modified_classes=len(mod),affected_support_pairs=count,critical_support_pairs=critical,distinct_restricted_sign_pairs=signpairs,max_critical_dot=maxdot,certificate_sha256=hashlib.sha256(SRC.read_bytes()).hexdigest(),scope='Independent brute-force65536-word kernels and direct critical sign-pair products; baseline/transfer replay separate.')
(HERE/'independent_verification.json').write_text(json.dumps(report,indent=2));print(report)
