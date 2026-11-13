"""Independent explicit sign-projection audit. No discovery/verifier imports."""
from pathlib import Path
import json,hashlib
import numpy as np
H=Path(__file__).resolve().parent; R=H.parents[3]; P=R/'verifications/improved/dim61-67-compatible-sign-families'
C=H.parent/'certificate.json'; B=P/'data/global64.json'; S=R/'verifications/improved/dim61-68-95-96-code-layers/source/dim61_96/code64_supports.npz'
c=json.loads(C.read_bytes());b=json.loads(B.read_bytes())
with np.load(S,allow_pickle=False) as z: parents=list(map(int,z['words']))
assert len(parents)==len(set(parents))==30828
assert all(0<s<2**64 and s.bit_count()==16 for s in parents)
old={r['support']:r for r in b['old']}; fixed={r['support']:r for r in b['new']}; modified={r['support']:r for r in c['modified_parent_classes']}; new={r['support']:r for r in c['new_classes']}
assert len(modified)==len(c['modified_parent_classes']) and len(new)==len(c['new_classes'])
assert set(modified)<=set(parents) and not set(modified)&set(old)
assert not set(new)&(set(parents)|set(fixed))
assert len(fixed)==len(b['new'])==1940
allrows=old|fixed|modified|new
X=np.arange(65536,dtype=np.uint32); pc=np.array([i.bit_count() for i in range(65536)],dtype=np.uint8)
cache={}; minweights={}
def signs(s):
 if s in cache:return cache[s]
 r=allrows[s]; assert 0<s<2**64 and s.bit_count()==16
 hs=r['checks']; phase=r['shift']; assert 0<=phase<65536
 assert all(type(h)==int and 0<=h<65536 for h in hs)
 ok=np.ones(65536,dtype=bool)
 for h in hs:ok &= (pc[X & h]%2==0)
 lin=X[ok];assert lin[0]==0 and len(lin)>1
 assert min(pc[lin[1:]])>=4
 assert 65535 in lin
 if s not in new:assert len(lin)==2048
 else:assert len(lin)==r['count']==2**r['dimension']
 vals=lin^phase;assert len(vals)==len(set(map(int,vals)))
 cache[s]=vals;minweights[str(s)]=int(min(pc[lin[1:]]));return vals
for s in modified:signs(s)
for s in new:signs(s)
# Modified old supports are safe against unchanged support classes independently of signs.
modpairs=0
for s in modified:
 for t in parents:
  if s!=t:assert (s&t).bit_count()<=8;modpairs+=1
 for t in fixed:assert (s&t).bit_count()<=8;modpairs+=1
# For high intersections, enumerate BOTH projected sign sets explicitly.
# At intersection9/10 any differing sign ensures dot=4*(m-2d)<=32.
critical=[];pairs=0
proj_cache={}
def project(s,inter):
 key=(s,inter)
 if key in proj_cache:return proj_cache[key]
 pos=[i for i in range(64) if s>>i&1]
 mask=sum(1<<j for j,i in enumerate(pos) if inter>>i&1)
 idx=[j for j in range(16) if mask>>j&1]
 vals=signs(s);out=np.zeros(len(vals),dtype=np.uint16)
 for k,j in enumerate(idx):out|=((vals>>j)&1).astype(np.uint16)<<k
 result=set(map(int,out));proj_cache[key]=result;return result
for i,s in enumerate(new):
 for t in parents+list(fixed)+list(new)[:i]:
  pairs+=1;inter=s&t;m=inter.bit_count()
  if m<=8:continue
  assert m<=10 and t in allrows
  A=project(s,inter);D=project(t,inter)
  assert A.isdisjoint(D),(s,t,m)
  critical.append({'left':s,'right':t,'intersection':m,'left_projection_size':len(A),'right_projection_size':len(D)})
# All sparse points have norm64, within-class dot<=64-8*4=32.
assert 16*2**2==64 and 64-8*4==32
# All non-sparse unaffected layers: dense +/-1, quadruple +/-4, axes +/-8.
assert 16*2<=32 and 4*4*2<=32 and 8*2<=32
# Owner classes contain only dense and weight4 points, no changed sparse point.
# Cross-cap product <=sqrt3/2*32=16sqrt3<32; poles vs sparse =0.
assert 16**2*3<32**2
base=2**28+(30828+1940)*2048+10416*16+128
assert c['base_count']==base==335711104
gain=sum(len(signs(s)) for s in new)
assert base+gain==c['total']==335749696
report={'status':'PASS','total64':base+gain,'total65':base+gain+68224+2,'gain':gain,'new_classes':len(new),'modified_classes':len(modified),'changed_old_support_pairs':modpairs,'new_support_pairs':pairs,'critical_pairs':critical,'explicit_sign_classes_enumerated':len(cache),'explicit_sign_words_enumerated':sum(map(len,cache.values())),'minimum_sign_weights':minweights,'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [C,B,S]},'scope':'New/modified classes independently enumerated. High-overlap pairs verified using disjoint explicit projected sign sets, not parity witness algorithm. Unchanged published construction and owner lift replayed separately.'}
(H/'report.json').write_text(json.dumps(report,indent=2)+'\n')
print({k:v for k,v in report.items() if k not in ['critical_pairs','minimum_sign_weights','sha256']})
