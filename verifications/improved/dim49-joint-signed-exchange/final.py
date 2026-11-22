"""Independent exact reconstruction and full final-pair check for tau(49)."""
import argparse,gzip,hashlib,json
from pathlib import Path
import numpy as np
DIR=Path(__file__).resolve().parent
SOURCE='ee848a6940bed96246dd3cdfa9e64eb10abd768757b40f86d7d6f114b1f54480'
TRADE='a88a834814f31d844f5ef095d9164640fa201bdf13d289a7fa8c4b553a12ed91'
def sha(b):return hashlib.sha256(b).hexdigest()
def verify(path):
 raw=path.read_bytes();assert sha(raw)==SOURCE
 rawdelta=(DIR/'trade.json.gz').read_bytes();assert sha(rawdelta)==TRADE
 j=json.loads(gzip.decompress(rawdelta))
 assert set(j)=={'removed_indices','added'} and all(type(v) is int for v in j['removed_indices'])
 assert all(len(v)==48 and all(type(a) is int and abs(a)<2**30 for a in v)for v in j['added'])
 R=list(j['removed_indices']);M=np.array(j['added'],'i8')
 a=np.load(path,allow_pickle=False);assert set(a.files)=={'points','gram'}
 S=a['points'];G=a['gram']
 assert S.dtype.kind in 'iu' and G.dtype.kind in 'iu'
 S=S.astype('i8');G=G.astype('i8')
 assert S.shape==(15622,48) and G.shape==(48,48) and np.array_equal(G,G.T)
 assert len(R)==len(set(R))==309 and all(0<=i<len(S)for i in R)
 assert M.shape==(343,48)
 assert not(set(map(tuple,M))&set(map(tuple,S)))
 present=np.ones(len(S),bool);present[R]=False
 C=np.r_[S[present],M];N=len(C)
 assert len({tuple(v)for v in C})==N==15656
 normmax=max(abs(int(C.min())),abs(int(C.max())))
 grammax=max(abs(int(G.min())),abs(int(G.max())))
 assert 48*normmax*grammax<2**63
 CG=C@G
 maxcg=max(abs(int(CG.min())),abs(int(CG.max())))
 assert 48*maxcg*normmax < 2**53
 assert np.all(np.sum(CG*C,axis=1)==6)
 Ct=C.T.astype(float)
 pairs=0
 for k in range(0,N,192):
  Z=CG[k:k+192].astype(float)@Ct
  for i,z in enumerate(Z):
   a=k+i
   assert z[a]==6 and np.all(z[:a]<=2)
   pairs+=a
 assert pairs==N*(N-1)//2
 return dict(status='PASS',dimension=49,baseline=52431624,added=len(M),removed=len(R),gain=N-len(S),selected=N,unordered_pairs=pairs,lower_bound=52416000+N+2,source_sha256=SOURCE,trade_sha256=TRADE)
def results():
    row = json.loads((DIR/'results.json').read_text())
    return [dict(dimension=row['dimension'], total=row['lower_bound'])]

if __name__=='__main__':
 if not __debug__:raise RuntimeError('Enable assertions')
 p=argparse.ArgumentParser();p.add_argument('--source',type=Path,default=DIR.parent/'dim49-60-p48-signed-refill/p48_signed_class_lns.npz')
 print(json.dumps(verify(p.parse_args().source),indent=2))
