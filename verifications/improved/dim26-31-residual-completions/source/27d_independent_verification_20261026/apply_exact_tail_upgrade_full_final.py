'''Apply a fully checked pure-tail code to a separately verified frozen nonpure base.'''
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import sys,json,hashlib,resource
from pathlib import Path
resource.setrlimit(resource.RLIMIT_AS,(256*1024**2,256*1024**2))
import numpy as np
D=Path(__file__).resolve().parent;P=D.parent/'27d_local_replacements_20261026'
which=int(sys.argv[1]);assert which in (30,31)
configs={30:(325,238796,72,3,'a0bc0ef49d0a163e84c7515d9b6787e28bec7675d1ca2af295c1cab92090138d','4c48503313b62b3aa5213afdfdc1bc3ac0fc5ad3522c5bd6e37468a0e5e0cb8d','14f522d79042b2672b1e445ada58573ed8390d51b3f395368e2a1bf3ccccbd4a'),31:(385,246596,126,2,'3da34aa68f3cb3ce99601b766fc0a170757b83651b025d81527452464c4e9a9f','ee86eaaa21092d4e686f4969bc4366eb01dcda79d81684255e53763017894c59','3a6760930e99381a2ece36fea83f78011d5b5ee43b4c5f587cfc5a32f9489f61')}
groups,total,count,rad,bs,es,ts=configs[which];m=which-24
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
bp=P/f'base{which}_groups{groups}.npz';ep=P/f'extension{which}_groups{groups}.npz';tp=P/f'pure_tail{which}_exceptional.npz'
assert [sha(bp),sha(ep),sha(tp)]==[bs,es,ts]
rp=D/f'verification{which}_{total}_full_replay.json';r=json.loads(rp.read_text())
assert r['status']=='PASS' and r['total']==total and r['base_sha256']==bs and r['extension_sha256']==es
old=np.load(bp,allow_pickle=False)['new_points'].astype(np.int32);free=np.load(ep,allow_pickle=False)['free_heads'].astype(np.int32)
assert old.shape==(groups*256,which) and free.shape==(m*96,which)
assert np.all(np.count_nonzero(old[:,24:],axis=1)==2) and np.all(np.isin(old[:,24:],[-2,0,2]))
assert np.all(np.count_nonzero(free[:,24:],axis=1)==1) and np.all(np.isin(free[:,24:],[-4,0,4]))
assert np.all(np.sum(old[:,:24]**2,axis=1)==24) and np.all(np.sum(free[:,:24]**2,axis=1)==16)
z=np.load(tp,allow_pickle=False);assert int(z['radicand'])==rad
assert z['rational'].shape==z['radical'].shape==(count,m)
C=[tuple(zip(map(int,a),map(int,b))) for a,b in zip(z['rational'],z['radical'])];assert len(set(C))==count

def nonnegative(a,b):
 if b==0:return a>=0
 if b>0:return a>=0 or rad*b*b>=a*a
 return a>=0 and a*a>=rad*b*b

def dot(x,y):return sum(a*c+rad*b*d for (a,b),(c,d) in zip(x,y)),sum(a*d+b*c for (a,b),(c,d) in zip(x,y))
for i,x in enumerate(C):
 assert dot(x,x)==(32,0)
 assert all(nonnegative(4-a,-b) and nonnegative(4+a,b) for a,b in x)
 for y in C[:i]:
  a,b=dot(x,y);assert nonnegative(16-a,-b)
answer=total-2*m*(m-1)+count
out={'status':'PASS_EXACT_UPGRADE','dimension':which,'verified_integer_base_count':total,'replacement_total':answer,'old_pure_count':2*m*(m-1),'new_pure_count':count,'base_sha256':bs,'extension_sha256':es,'pure_tail_sha256':ts,'integer_verification_report_sha256':sha(rp),'old_tail_pattern':'exactly two±2','free_tail_pattern':'exactly one±4','nonpure_head_norms_squared':[24,16],'proof':'cube bound4 implies all new/old and new/free products≤16; equator products0; exact pure pairs/norms checked; nonzero heads exclude duplicates','scope':'Depends on fresh full integer verification from the accompanying standalone checker and matching frozen inputs; no earlier PASS report or compressed compatibility is inherited.'}
(D/f'verification{which}_{answer}_exceptional_full_replay.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
