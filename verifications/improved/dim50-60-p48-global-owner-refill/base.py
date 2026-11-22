from pathlib import Path
import numpy as np,gzip,json,hashlib
R=Path('/workspace/scratch/941a2e823e69/research')
T=R/'p48_pole_split_20261122/combined'
F=R/'p48_signed_flow_compact_20261122'
B0=R/'mid_49_69_20261121/verifications/improved'
POLE=[6,12,24,40,72,126,240,306,510,604,841]
PREV=[52462872,52509740,52603402,52728054,52976623,53393546,54263398,54760693,56270076,56952414,58631571]
TH={
54:('d54','witness.json','8b62dc85b2d16c339d186d899b1179021a16621358f4b270d07bdb0352c4f292','head_triple',156),
55:('d55trine','witness.json','ff2d61cdba830024e54a9ab371ddc63ff256fa27a7b68d5c4323925edc733fd0','head_triples',126),
57:('d57trine','witnesses.json','08c0cb6df122d888258b7e0d050766c97809b226980181dcda2516f977386b5c','head_triples',264),
59:('d59trine','witnesses.json','00ba9207ec2f3a39616667255c645708d9798dca408114b0b9eb6a8233bf94a5','head_triples',168)}
def SHA(s):return hashlib.sha256(s).hexdigest()
def m(a):return max(abs(int(a.min())),abs(int(a.max()))) if a.size else 0
def fp(A,B):
 assert A.ndim==B.ndim==2 and A.shape[1]==B.shape[0] and A.shape[1]*m(A)*m(B)<2**53
 return A.astype('d')@B.astype('d')
def ip(A,B):return fp(A,B).astype('i8')
def key(A):
 assert A.ndim==2 and A.shape[1]==48 and m(A)<2**15
 return np.ascontiguousarray(A,dtype='<i2').view('V96').ravel()
def matched(sorted_keys,queries):
 q=np.searchsorted(sorted_keys,queries)
 keep=q<len(sorted_keys);ks=np.flatnonzero(keep);keep[ks]=sorted_keys[q[ks]]==queries[ks]
 return q,keep
def tails(d,TP=T):
 if d==60:
  raw=(TP/'d60cube/witness.json').read_bytes()
  assert SHA(raw)=='15316f8280336477cea06d80fcff9132af39ca09313ceebb37657ba1eeae78b1'
  H=np.array(json.loads(raw),'i8')
  return np.r_[H,-H],3583
 if d not in TH:return np.empty((0,48),'i8'),0
 u,v,sha,tag,gain=TH[d];s=(TP/u/v).read_bytes();assert SHA(s)==sha
 x=json.loads(s)
 if isinstance(x,list):x=next(w for w in x if w['dimension']==d)
 H=np.array(x[tag],'i8').reshape(-1,48)
 return np.r_[H,-H],gain
class Problem:
 def __init__(self,d,root=B0):
  bp=root/'dim49-60-p48-pointwise-refill';sp=root/'dim49-60-p48-signed-refill'
  p=bp/f'heads{d}.npz';q=sp/f'p48_signed_delta{d}.json.gz'
  raw=q.read_bytes();D=json.loads(gzip.decompress(raw))
  self.ph=SHA(p.read_bytes());self.oh=SHA(raw)
  assert D['bound']==PREV[d-50] and D['dimension']==d
  V=np.load(p,allow_pickle=False)
  self.g=V['gram'].astype('i8');self.T=V['transforms'].astype('i8')
  self.base=V['base'].astype('i8');self.bit=np.unpackbits(V['retained_bits'],axis=1)[:,:len(self.base)].astype(bool)
  self.X=V['extra_points'].astype('i8');self.ix=V['extra_offsets']
  self.trade={v['class_index']:v for v in D['changes']};self.d=d
 def rows(self,i):
  X=np.r_[self.base[self.bit[i]],self.X[self.ix[i]:self.ix[i+1]]]
  Y=np.r_[X,-X]
  if i in self.trade:
   t=self.trade[i];drop=np.sort(key(np.array(t['remove'],'i8').reshape(-1,48)));g=key(Y)
   if len(drop):
    _,R=matched(drop,g);assert int(R.sum())==len(drop)
    Y=Y[~R]
   Y=np.r_[Y,np.array(t['add'],'i8').reshape(-1,48)]
  assert len(Y)<32768
  return Y
