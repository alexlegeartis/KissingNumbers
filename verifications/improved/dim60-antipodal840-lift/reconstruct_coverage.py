"""Rebuild deterministic full-image maximum coverage (not needed for verification)."""
from pathlib import Path
import numpy as np,json
D=Path(__file__).resolve().parent;z=np.load(D/'source/family_refilled.npz');base=z['base'].astype(np.int64);G=z['gram'].astype(np.int64);hist=np.load(D/'source/family_55_61.npz')['transforms'].astype(np.int64);gens=np.load(D/'source/p48p_gens.npy').astype(np.int64)
T=[];X=[]
for i,ref in enumerate(z['transform_references']):
 if ref>=0:t=hist[ref]
 else:
  t=np.eye(48,dtype=np.int64)
  for k in z['transform_words'][i]:t=t@gens[k]
 assert np.array_equal(t@G@t.T,G);T.append(t);x=base@t;j=np.argmax(x!=0,axis=1);x*=np.sign(x[np.arange(len(x)),j])[:,None];assert np.max(abs(x))<128;X.append(x.astype(np.int8))
X=np.concatenate(X);keys=X.view(np.dtype((np.void,48))).ravel();uniq,inv=np.unique(keys,return_inverse=True);ids=inv.reshape(577,7783);del X,keys,inv
transforms=np.array(T)
chosen=np.load(D/'selected_heads_full.npz')['class_ids'].copy();M=int(ids.max())+1;history=[]
while True:
 counts=np.bincount(ids[chosen].ravel(),minlength=M);owner=np.full(M,-1,dtype=np.int32)
 for j,i in enumerate(chosen):owner[ids[i][counts[ids[i]]==1]]=j
 loss=np.array([np.sum(counts[ids[i]]==1) for i in chosen]);best=(0,-1,-1)
 for add in sorted(set(range(577))-set(chosen)):
  vv=ids[add];new=int(np.sum(counts[vv]==0));oo=owner[vv];recover=np.bincount(oo[oo>=0],minlength=420);net=new-loss+recover;drop=int(np.argmax(net));gain=int(net[drop])
  if gain>best[0]:best=(gain,add,drop)
 if best[0]<=0:break
 gain,add,drop=best;old=int(chosen[drop]);chosen[drop]=add;history.append(dict(gain=gain,add=add,remove=old));print(len(history),history[-1],flush=True)
 if len(history)>=200:break
_,first=np.unique(ids[chosen].ravel(),return_index=True);bits=np.zeros(420*7783,np.uint8);bits[first]=1;bits=bits.reshape(420,7783);S=len(first)
np.savez_compressed(D/'optimized420.npz',class_ids=chosen,transforms=transforms[chosen],retained_bits=np.packbits(bits,axis=1),sizes=bits.sum(1))
r=dict(selected_lines=S,bound60=52416000+2*S+841,gain_over_published=2*(S-3062990),exchanges=history,one_exchange_local_optimum=best[0]<=0)
(D/'coverage_exchanges.json').write_text(json.dumps(r,indent=2));print(r,flush=True)
