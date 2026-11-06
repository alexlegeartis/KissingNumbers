from pathlib import Path
import numpy as np,json,hashlib
if not __debug__:raise SystemExit('Run without -O')
def mul(a,b):
 assert a.shape[1]*int(np.max(abs(a)))*int(np.max(abs(b)))<2**63
 return a@b
D=Path(__file__).resolve().parent;S=D/'source';z=np.load(S/'family_refilled.npz');base=z['base'].astype(np.int64);G=z['gram'].astype(np.int64);a=np.load(D/'optimized420.npz');bits=np.unpackbits(a['retained_bits'],axis=1);assert not bits[:,7783:].any();bits=bits[:,:7783];assert np.array_equal(bits.sum(1),a['sizes']);assert len(bits)==420
assert G.shape==(48,48) and np.array_equal(G,np.load(S/'p48p_gram.npy'));assert np.all(np.sum(mul(base,G)*base,axis=1)==6)
for start in range(0,7783,256):
 h=mul(mul(base[start:start+256],G),base.T)
 for j,row in enumerate(h):
  i=start+j;assert np.all(abs(row[:i])<=2) and np.all(abs(row[i+1:])<=2)
X=[]
for t,mask in zip(a['transforms'].astype(np.int64),bits):
 assert np.array_equal(mul(mul(t,G),t.T),G);x=mul(base[mask.astype(bool)],t);assert np.all(np.sum(mul(x,G)*x,axis=1)==6);j=np.argmax(x!=0,axis=1);x*=np.sign(x[np.arange(len(x)),j])[:,None];assert np.max(abs(x))<128;X.append(x.astype(np.int8))
X=np.concatenate(X);keys=X.view(np.dtype((np.void,48))).ravel();assert len(np.unique(keys))==len(keys)==3063342
r=dict(status='PASS',selected_lines=len(keys),classes=420,mother_lines=len(keys),bound60=52416000+2*len(keys)+841,gain_over_58542821=2*(len(keys)-3062990),base_pairs_checked=7783*7782//2,verified_class_isometries=420,verification='Full exact coordinate disjointness; tail/pole and mother-lattice premises from published f0083954 package',sha256=hashlib.sha256((D/'optimized420.npz').read_bytes()).hexdigest())
(D/'coverage_verification.json').write_text(json.dumps(r,indent=2));print(r)
