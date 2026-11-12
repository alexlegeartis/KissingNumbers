"""Read-only complete signed-class check. --cpp adds independent scalar int64 replay."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,json,struct,subprocess,tempfile,time
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
PARENT=ROOT/'verifications/improved/dim49-60-p48-pointwise-refill'
GRAM_SHA='01a5d107e4424092e5300e201592aaff62804ee27eeb1d23b6e516249897d3e1'

def check(path,cpp=False):
 raw=Path(path).read_bytes();a=np.load(path,allow_pickle=False)
 assert set(a.files)=={'points','gram'}
 X=a['points'];G=a['gram'];assert X.dtype.kind in 'iu' and G.dtype.kind in 'iu'
 X=X.astype(np.int64);G=G.astype(np.int64)
 assert X.ndim==2 and X.shape[1]==48 and G.shape==(48,48)
 gram_path=PARENT/'source/p48p_gram.npy';assert hashlib.sha256(gram_path.read_bytes()).hexdigest()==GRAM_SHA
 assert np.array_equal(G,np.load(gram_path))
 assert len(np.unique(X,axis=0))==len(X)
 mx=int(np.max(np.abs(X)));mg=int(np.max(np.abs(G)))
 assert 48*mx*mg<2**63
 XG=X@G;md=int(np.max(np.abs(XG)));assert 48*mx*md<2**53
 assert np.all(np.sum(XG*X,axis=1)==6)
 pairs=0;maximum=-6
 for i in range(0,len(X),256):
  dots=XG[i:i+256].astype(np.float64)@X.astype(np.float64).T
  for j,row in enumerate(dots):
   at=i+j;assert row[at]==6
   if at:maximum=max(maximum,int(row[:at].max()));assert np.all(row[:at]<=2)
   pairs+=at
 assert pairs==len(X)*(len(X)-1)//2
 keys=set(map(tuple,X));asym=sum(tuple(-x) not in keys for x in X)
 # The proof of the 49D lift works for signed heads separately.
 assert Q(3,4)+Q(1,4)==1
 assert Q(3,4)*Q(1,3)+Q(1,4)==Q(1,2)
 assert Q(3,4)-Q(1,4)==Q(1,2)
 assert Q(3,4)*Q(1,2)**2<Q(1,2)**2
 result=dict(status='PASS',signed_points=len(X),unordered_pairs=pairs,maximum_product=maximum,points_without_antipodes=asym,kissing49=52416000+len(X)+2,gain49=len(X)-15598,certificate_sha256=hashlib.sha256(raw).hexdigest(),published_gram_sha256=GRAM_SHA,arithmetic='integer matrix products followed by binary64 integer products; every product and partial sum strictly below2^53',inherited_inputs=['The pinned P48p Gram is an even unimodular48D lattice of minimum6.','Its complete minimal shell has52,416,000 vectors.'],tail='Two poles in one dimension; each chosen signed mother point is replaced by its two opposite-tail lifts.')
 if cpp:
  with tempfile.TemporaryDirectory() as t:
   t=Path(t);data=t/'vectors.bin';data.write_bytes(struct.pack('<Q',len(X))+G.astype('<i8').tobytes()+X.astype('<i8').tobytes());exe=t/'check'
   subprocess.run(['g++','-O3','-std=c++17',str(HERE/'check_signed_class.cpp'),'-o',str(exe)],check=True)
   out=subprocess.run([str(exe),str(data)],check=True,capture_output=True,text=True)
   result['independent_cpp']=json.loads(out.stdout)
   assert result['independent_cpp']['signed_points']==len(X) and result['independent_cpp']['unordered_pairs']==pairs
 return result
if __name__=='__main__':
 if not __debug__:raise SystemExit('Assertions must be enabled.')
 parser=argparse.ArgumentParser();parser.add_argument('certificate',nargs='?',default=str(HERE/'p48_signed_class_lns.npz'));parser.add_argument('--cpp',action='store_true');a=parser.parse_args();print(json.dumps(check(a.certificate,a.cpp),indent=2))
