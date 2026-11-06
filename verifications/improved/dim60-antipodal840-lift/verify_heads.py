"""Exact replay of 420 generator-word P48 classes, without historical data."""
import numpy as np,json,hashlib,time,sys
FULL='--full' in sys.argv
from pathlib import Path
if not __debug__:raise SystemExit('Run without -O: exact verification uses assertions.')
D=Path(__file__).parent;S=D/'source';start=time.monotonic()
z=np.load(S/'family_refilled.npz');base=z['base'].astype(np.int64);G=z['gram'].astype(np.int64);gens=np.load(S/'p48p_gens.npy').astype(np.int64)
def mul(a,b):
 assert a.shape[1]*int(np.max(abs(a)))*int(np.max(abs(b)))<2**63
 return a@b
def canon(X):
 assert np.all(np.any(X!=0,axis=1));j=np.argmax(X!=0,axis=1);return X*np.sign(X[np.arange(len(X)),j])[:,None]
assert G.shape==(48,48) and np.array_equal(G,G.T)
assert np.array_equal(G,np.load(S/'p48p_gram.npy'))
assert base.shape==(7783,48) and len(np.unique(canon(base),axis=0))==7783
BG=mul(base,G);assert np.all((BG*base).sum(1)==6)
checked=0
for st in range(0,len(base),256):
 H=mul(BG[st:st+256],base.T)
 for k,row in enumerate(H):
  i=st+k;assert np.all(abs(row[:i])<=2) and np.all(abs(row[i+1:])<=2);checked+=i
assert checked==30283653
print('base pairs pass',checked,flush=True)
for a in gens:assert np.array_equal(mul(mul(a,G),a.T),G)
refs=z['transform_references'];ids=np.arange(577) if FULL else np.flatnonzero(refs<0);ids=ids[np.argsort(-z['sizes'][ids],kind='stable')][:420];assert len(ids)==420
historical=np.load(S/'family_55_61.npz')['transforms'].astype(np.int64) if FULL else None
bits=np.unpackbits(z['retained_bits'],axis=1);assert not bits[:,len(base):].any();bits=bits[:,:len(base)];assert np.array_equal(bits.sum(1),z['sizes'])
coords=[];trans=[];sizes=[]
for count,i in enumerate(ids):
 word=z['transform_words'][i]
 if refs[i]>=0:
  assert FULL and np.all(word==-1);T=historical[refs[i]]
 else:
  assert refs[i]==-1 and np.all((word>=0)&(word<len(gens)))
  T=np.eye(48,dtype=np.int64)
  for k in word:T=mul(T,gens[k])
 assert np.array_equal(mul(mul(T,G),T.T),G)
 X=canon(mul(base[bits[i].astype(bool)],T));assert np.max(abs(X))<128
 assert np.all((mul(X,G)*X).sum(1)==6)
 coords.append(X.astype(np.int8));trans.append(T);sizes.append(len(X))
 if count%50==0:print('classes',count+1,'seconds',time.monotonic()-start,flush=True)
X=np.concatenate(coords);keys=np.ascontiguousarray(X).view(np.dtype((np.void,48))).ravel();assert len(np.unique(keys))==len(keys)
selected_lines=len(X);assert selected_lines==(3062990 if FULL else 3035307)
report={'status':'PASS','base_lines':len(base),'base_pairs_checked':checked,'classes':420,'selected_classes':ids.tolist(),'selected_lines':selected_lines,'all_selected_transforms_isometric':True,'line_disjointness':'full canonical 48-coordinate byte comparison','mother_shell_size':52416000,'mother_minimum_squared_norm':6,'mother_shell_scope':'Inherited published P48p minimum and shell population; not re-enumerated.','tail_antipodal_pairs':420,'pole_count':841,'bound60':52416000+2*selected_lines+841,'comparison_bound60':57968343,'gain':2*(selected_lines-2775751),'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [S/'family_refilled.npz',S/'p48p_gens.npy']}}
np.savez_compressed(D/('selected_heads_full.npz' if FULL else 'selected_heads.npz'),class_ids=ids,transforms=np.array(trans),sizes=np.array(sizes))
(D/('head_verification_full.json' if FULL else 'head_verification.json')).write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k!='selected_classes'}),flush=True)
