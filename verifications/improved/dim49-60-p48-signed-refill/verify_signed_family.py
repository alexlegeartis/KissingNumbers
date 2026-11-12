"""Read-only full P48 signed-family geometry and global ownership replay.
The existing base's entire absolute-product matrix is checked once. For each
final class, every exceptional/retained-base and exceptional/exceptional pair
is checked, allowing individual signs of base vectors to be removed.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,gzip,hashlib,importlib.util,json,time
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
PARENT=ROOT/'verifications/improved/dim49-60-p48-pointwise-refill'
GRAM_SHA='01a5d107e4424092e5300e201592aaff62804ee27eeb1d23b6e516249897d3e1'
POLES=[2,6,12,24,40,72,126,240,306,510,604,841]
CLASSES=[1,3,6,12,20,36,63,120,153,255,302,420]

def integer(a):
 assert a.dtype.kind in 'iu';assert not a.size or max(abs(int(a.min())),abs(int(a.max())))<2**30
 return a.astype(np.int64)
def mag(a):return max(abs(int(a.min())),abs(int(a.max()))) if a.size else 0
def iproduct(a,b):
 assert a.ndim==b.ndim==2 and a.shape[1]==b.shape[0]
 assert a.shape[1]*mag(a)*mag(b)<2**63
 return a@b
def fproduct(a,b):
 assert a.ndim==b.ndim==2 and a.shape[1]==b.shape[0]
 assert a.shape[1]*mag(a)*mag(b)<2**53
 return a.astype(np.float64)@b.astype(np.float64)
def canonical(a):
 assert a.ndim==2 and a.shape[1]==48
 if not len(a):return a.copy(),np.empty(0,dtype=np.int64)
 assert np.all(np.any(a!=0,axis=1))
 jj=np.argmax(a!=0,axis=1);sign=np.sign(a[np.arange(len(a)),jj]);return a*sign[:,None],sign
def keys(a):
 assert mag(a)<2**15
 return np.ascontiguousarray(a,dtype=np.int16).view(np.dtype((np.void,96))).ravel()
def signed_keys(a):return [row.astype('<i2').tobytes() for row in a]
def check_norm(a,G):
 q=iproduct(a,G);assert a.shape[1]*mag(q)*mag(a)<2**63
 assert np.all(np.sum(q*a,axis=1)==6)
 return q

def check_base():
 path=PARENT/'source/p48p_gram.npy';assert hashlib.sha256(path.read_bytes()).hexdigest()==GRAM_SHA
 G=integer(np.load(path,allow_pickle=False));assert G.shape==(48,48) and np.array_equal(G,G.T)
 f=np.load(PARENT/'heads49.npz',allow_pickle=False);B=integer(f['base']);assert B.shape==(7783,48)
 BG=check_norm(B,G);bc,_=canonical(B);assert len(np.unique(keys(bc)))==len(B)
 count=0
 for st in range(0,len(B),256):
  H=fproduct(BG[st:st+256],B.T)
  for j,h in enumerate(H):
   i=st+j;assert h[i]==6;assert np.all(abs(h[:i])<=2) and np.all(abs(h[i+1:])<=2);count+=i
 return B,G,count

def verify_dimension(d,B,G):
 start=time.time();path=HERE/f'p48_signed_delta{d}.json.gz';raw=path.read_bytes();decoded=gzip.decompress(raw);D=json.loads(decoded)
 source=PARENT/f'heads{d}.npz';a=np.load(source,allow_pickle=False)
 assert np.array_equal(a['base'],B) and np.array_equal(a['gram'],G)
 T=integer(a['transforms']);m=len(T);assert m==CLASSES[d-49] and T.shape==(m,48,48)
 bits=a['retained_bits'];assert bits.dtype==np.uint8 and bits.shape==(m,973)
 bits=np.unpackbits(bits,axis=1);assert not bits[:,len(B):].any();bits=bits[:,:len(B)].astype(bool)
 E=integer(a['extra_points']);O=integer(a['extra_offsets']);assert O.shape==(m+1,) and O[0]==0 and O[-1]==len(E) and np.all(np.diff(O)>=0)
 assert E.ndim==2 and E.shape[1]==48;check_norm(E,G)
 changes={int(c['class_index']):c for c in D['changes']};assert len(changes)==len(D['changes'])
 assert all(0<=i<m for i in changes)
 Blookup={tuple(x):(j,1) for j,x in enumerate(B)};Blookup.update({tuple(-x):(j,-1) for j,x in enumerate(B)})
 original_globals=[];deleted_globals=[];added_globals=[];class_sizes=[];checked=0;removed_total=0;added_total=0
 for i in range(m):
  t=T[i];assert np.array_equal(iproduct(iproduct(t,G),t.T),G)
  oldextra=E[O[i]:O[i+1]];oldlines=np.r_[B[bits[i]],oldextra]
  oldkeys={tuple(x) for x in oldlines}|{tuple(-x) for x in oldlines}
  assert len(oldkeys)==2*len(oldlines)
  pos=bits[i].copy();neg=bits[i].copy();exceptions=np.r_[oldextra,-oldextra]
  if i in changes:
   c=changes[i];remove=integer(np.array(c['remove'],dtype=np.int64).reshape(-1,48));add=integer(np.array(c['add'],dtype=np.int64).reshape(-1,48))
   remkeys=set(map(tuple,remove));addkeys=set(map(tuple,add))
   assert len(remkeys)==len(remove) and len(addkeys)==len(add)
   assert remkeys<=oldkeys and not(addkeys & (oldkeys-remkeys))
   for x in remove:
    pair=Blookup.get(tuple(x))
    if pair:
     j,s=pair
     if s>0:assert pos[j];pos[j]=False
     else:assert neg[j];neg[j]=False
   exceptions=np.array([x for x in exceptions if tuple(x) not in remkeys],dtype=np.int64).reshape(-1,48)
   exceptions=np.r_[exceptions,add]
   deleted_globals.append(iproduct(remove,t));added_globals.append(iproduct(add,t))
   removed_total+=len(remove);added_total+=len(add)
  # All signed base/base pairs follow from the freshly checked absolute bound.
  EG=check_norm(exceptions,G)
  for st in range(0,len(exceptions),256):
   H=fproduct(EG[st:st+256],B.T)
   assert np.all(H[:,pos]<=2) and np.all(H[:,neg]>=-2)
   checked+=len(H)*(int(pos.sum())+int(neg.sum()))
  for st in range(0,len(exceptions),256):
   H=fproduct(EG[st:st+256],exceptions.T)
   for j,h in enumerate(H):
    at=st+j;assert h[at]==6;assert np.all(h[:at]<=2) and np.all(h[at+1:]<=2);checked+=at
  global_lines=iproduct(oldlines,t);canon,_=canonical(global_lines)
  original_globals.append(canon.astype(np.int16))
  class_sizes.append(int(pos.sum())+int(neg.sum())+len(exceptions))
 original=np.concatenate(original_globals);del original_globals
 source_keys=keys(original);sortedkeys=np.sort(source_keys)
 # Complete parent-coordinate ownership; numerical hashes are not used.
 assert not np.any(sortedkeys[1:]==sortedkeys[:-1])
 removed=np.concatenate(deleted_globals) if deleted_globals else np.empty((0,48),dtype=np.int64)
 added=np.concatenate(added_globals) if added_globals else np.empty((0,48),dtype=np.int64)
 # Independent signed complete-coordinate records settle orientation exactly.
 deleted=set(signed_keys(removed));assert len(deleted)==len(removed)
 inserted=set(signed_keys(added));assert len(inserted)==len(added)
 for vectors,require_old in ((removed,True),(added,False)):
  canon,sign=canonical(vectors);query=keys(canon);locations=np.searchsorted(sortedkeys,query)
  for x,q,j in zip(vectors,query,locations):
   owned=j<len(sortedkeys) and sortedkeys[j]==q
   if require_old:assert owned
   elif owned:assert x.astype('<i2').tobytes() in deleted
 selected=2*len(original)-len(removed)+len(added)
 assert selected==sum(class_sizes)
 bound=52416000+selected+POLES[d-49]
 assert D['bound']==bound and D['baseline']==52416000+2*len(original)+POLES[d-49]
 assert D['gain']==len(added)-len(removed) and D['dimension']==d
 return dict(status='PASS',dimension=d,bound=bound,baseline=D['baseline'],gain=D['gain'],classes=m,changed_classes=len(changes),selected_signed_mother_points=selected,removed_signed_points=removed_total,added_signed_points=added_total,explicit_exceptional_pairs=checked,global_original_lines=len(original),global_added_signed_points=len(inserted),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),delta_sha256=hashlib.sha256(raw).hexdigest(),delta_decoded_sha256=hashlib.sha256(decoded).hexdigest(),seconds=time.time()-start)

if __name__=='__main__':
 if not __debug__:raise SystemExit('Assertions must be enabled.')
 parser=argparse.ArgumentParser();parser.add_argument('dimensions',nargs='*',type=int);parser.add_argument('--tails',action='store_true');args=parser.parse_args()
 B,G,pairs=check_base();reports=[]
 for d in args.dimensions or range(49,61):
  assert 49<=d<=60
  report=verify_dimension(d,B,G);reports.append(report);print(json.dumps(report),flush=True)
 summary=dict(status='PASS',base_pairs=pairs,dimensions=reports,arithmetic='Every int64 and binary64 multiplication has an explicit absolute partial-sum bound; no tolerance is used.',inherited_inputs=['Published P48p Gram identification, lattice minimum6 and52,416,000 minimal vectors.'],head_scope='Complete base pairs, all final exceptional/retained-base and exceptional/exceptional pairs, integral isometries and full-coordinate signed ownership are freshly checked.')
 if args.tails:
  spec=importlib.util.spec_from_file_location('tail_check',PARENT/'tail_and_lift.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
  summary['tails']=module.verify()
  summary['signed_lift_count']='M - S + 2S + P = 52416000 + S + P; classes need not be antipodal.'
 print(json.dumps({'summary':summary}),flush=True)
