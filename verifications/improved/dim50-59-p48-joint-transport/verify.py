import argparse,gzip,hashlib,json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
G_SHA='01a5d107e4424092e5300e201592aaff62804ee27eeb1d23b6e516249897d3e1'
POLES=[6,12,24,40,72,126,240,306,510,604,841]
DIMS=list(range(50,61))
PREV=[52462872,52509740,52603402,52728054,52976623,53393546,54263398,54760693,56270076,56952414,58631571]
def sha(b):return hashlib.sha256(b).hexdigest()
def mag(a):return max(abs(int(a.min())),abs(int(a.max())))if a.size else 0
def integer(a):
 a=np.asarray(a);assert a.dtype.kind in 'iu' and mag(a)<2**30;return a.astype('i8')
def product(a,b):
 assert a.ndim==b.ndim==2 and a.shape[1]==b.shape[0] and a.shape[1]*mag(a)*mag(b)<2**63
 return a@b
def fproduct(a,b):
 assert a.ndim==b.ndim==2 and a.shape[1]==b.shape[0] and a.shape[1]*mag(a)*mag(b)<2**53
 return a.astype(float)@b.astype(float)
def canonical(a):
 assert a.ndim==2 and a.shape[1]==48
 j=np.argmax(a!=0,axis=1);s=np.sign(a[np.arange(len(a)),j]);assert np.all(abs(s)==1)
 return a*s[:,None]
def rowkeys(a):
 assert mag(a)<2**15
 return np.ascontiguousarray(a,dtype='<i2').view('V96').ravel()
def rowkey(a):
 assert mag(a)<2**15
 return np.asarray(a,dtype='<i2').tobytes()
def checknorm(a,G):
 q=product(a,G);assert 48*mag(q)*mag(a)<2**63
 assert np.all(np.sum(q*a,axis=1)==6)
 return q
COMMON='a88a834814f31d844f5ef095d9164640fa201bdf13d289a7fa8c4b553a12ed91'
def verify(d,P,S,manifest):
 rawcommon=(HERE/'common_trade.json.gz').read_bytes();assert sha(rawcommon)==COMMON
 common=integer(json.loads(gzip.decompress(rawcommon))['added']);assert common.shape==(343,48)
 p=manifest[str(d)];a=HERE/f'transfer{d}.json.gz';r=a.read_bytes();assert sha(r)==p['transport_sha256']
 q=json.loads(gzip.decompress(r));assert q['dimension']==d and q['previous']==PREV[d-50]
 old=(S/f'p48_signed_delta{d}.json.gz').read_bytes()
 assert sha(old)==p['signed_delta_sha256']==q['prior_delta_sha256']
 z=json.loads(gzip.decompress(old));assert z['dimension']==d and z['bound']==q['previous']
 op=P/f'heads{d}.npz';assert sha(op.read_bytes())==p['pointwise_sha256']==q['original_sha256']
 gp=P/'source/p48p_gram.npy';assert sha(gp.read_bytes())==G_SHA
 F=np.load(op,allow_pickle=False);G=integer(F['gram']);B=integer(F['base']);T=integer(F['transforms'])
 assert G.shape==(48,48) and np.array_equal(G,np.load(gp,allow_pickle=False)) and B.shape==(7783,48)
 ori=np.unpackbits(F['retained_bits'],axis=1);assert not ori[:,len(B):].any();ori=ori[:,:len(B)].astype(bool)
 ext=integer(F['extra_points']);off=integer(F['extra_offsets'])
 assert len(ori)==len(T)==len(off)-1 and off[0]==0 and off[-1]==len(ext) and np.all(np.diff(off)>=0)
 oldd={int(c['class_index']):c for c in z['changes']};curd={int(c['class_index']):c for c in q['changes']}
 assert len(oldd)==len(z['changes']) and len(curd)==len(q['changes'])
 assert all(0<=c<len(T)for c in oldd) and all(0<=c<len(T)for c in curd)
 orig=[];oldextra=set();new=set();newcanon=[];removed=added=oldsize=0;pairs=0
 for i,t in enumerate(T):
  Y=np.r_[B[ori[i]],ext[off[i]:off[i+1]]]
  orig.append(rowkeys(canonical(product(Y,t))))
  oldsize+=2*len(Y)
  if i not in oldd and i not in curd:continue
  assert np.array_equal(product(product(t,G),t.T),G)
  Y=np.r_[Y,-Y]
  if i in oldd:
   c=oldd[i];rm=set(map(tuple,c['remove']));ob=np.asarray(c['add'],dtype='i8').reshape(-1,48)
   assert len(rm)==len(c['remove'])
   oldextra.update(rowkey(v)for v in product(ob,t))
   Y=np.r_[np.array([v for v in Y if tuple(v)not in rm],dtype='i8').reshape(-1,48),ob]
   oldsize+=len(ob)-len(rm)
  if i not in curd:continue
  c=curd[i];ix=c['insert'];assert all(type(j) is int and 0<=j<343 for j in ix) and len(set(ix))==len(ix)
  W=common[ix];A=checknorm(W,G)
  hits=fproduct(A,Y.T)>2;mask=hits.any(axis=0);rmcount=int(mask.sum())
  Y1=Y[~mask]
  D=fproduct(A,W.T);assert np.all(D.diagonal()==6) and np.all(D[np.triu_indices(len(W),1)]<=2)
  pairs+=len(W)*len(Y1)+len(W)*(len(W)-1)//2
  Gl=product(W,t)
  newcanon.append(rowkeys(canonical(Gl)))
  sk={rowkey(v)for v in Gl};assert len(sk)==len(W) and not(new&sk)
  new.update(sk)
  added+=len(W);removed+=rmcount
 K=np.sort(np.concatenate(orig));del orig
 assert np.all(K[1:]!=K[:-1])
 J=np.concatenate(newcanon)
 I=np.searchsorted(K,J);k=I<len(K)
 assert np.all(K[I[k]]!=J[k]) and not (new&oldextra)
 assert oldsize+52416000+POLES[d-50]==q['previous']
 assert added==q['added'] and removed==q['removed'] and added-removed==q['gain']>0
 assert q['bound']==q['previous']+q['gain']
 return dict(status='PASS',dimension=d,baseline=q['previous'],lower_bound=q['bound'],gain=q['gain'],added=added,removed=removed,
  changed_classes=len(curd),exact_new_pairs=pairs,transport_sha256=p['transport_sha256'])
if __name__=='__main__':
 if not __debug__:raise RuntimeError('Assertions required')
 p=argparse.ArgumentParser();p.add_argument('dimensions',type=int,nargs='*')
 p.add_argument('--pointwise',type=Path,default=HERE.parent/'dim49-60-p48-pointwise-refill')
 p.add_argument('--signed',type=Path,default=HERE.parent/'dim49-60-p48-signed-refill')
 a=p.parse_args()
 m=json.loads((HERE/'sources.json').read_bytes())
 for d in a.dimensions or sorted(map(int,m)):
  print(json.dumps(verify(d,a.pointwise,a.signed,m)),flush=True)
