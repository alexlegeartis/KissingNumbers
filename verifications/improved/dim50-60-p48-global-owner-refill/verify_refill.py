"""Exact new ownership / tail / two-stage P48 refill verification."""
from base import *
from pathlib import Path
import argparse,sys
HERE=Path(__file__).resolve().parent
ADV='953da203bd18bb961e7b558c0f63c27fb743ecb9a7c617eb133307d3c77ed935'
REF='ee848a6940bed96246dd3cdfa9e64eb10abd768757b40f86d7d6f114b1f54480'
GP='01a5d107e4424092e5300e201592aaff62804ee27eeb1d23b6e516249897d3e1'
def advanced():
 p=HERE/'advanced60_witness.json';z=p.read_bytes();assert SHA(z)==ADV
 a=json.loads(z)
 W=np.concatenate([np.array(a['frame'],'i8'),np.array(a['sparse'],'i8').reshape(-1,48),*[np.array(r['tetrahedra'],'i8').reshape(-1,48)for r in a['upgrade']]])
 # The pinned advanced tail witness independently proves its 4,799-point gain.
 assert W.shape==(64,48)
 return np.r_[W,-W],4799
def run(a, overrides=None, expected=None):
 overrides={} if overrides is None else overrides
 BF=(HERE/'bank.json.gz').read_bytes();A=np.array(json.loads(gzip.decompress(BF)),'i8')
 assert A.shape==(4000,48) and len(np.unique(key(A)))==len(A)
 refs=a.siblings/'dim49-60-p48-signed-refill/p48_signed_class_lns.npz';assert SHA(refs.read_bytes())==REF
 gp=a.siblings/'dim49-60-p48-pointwise-refill/source/p48p_gram.npy';assert SHA(gp.read_bytes())==GP
 with np.load(refs,allow_pickle=False)as x:
  Ref=x['points'].astype('i8');G=x['gram'].astype('i8')
 assert np.array_equal(G,np.load(gp,allow_pickle=False))
 n=len(A);r0=np.argsort(key(Ref));rk=key(Ref[r0]);assert len(np.unique(rk))==len(rk)
 M=json.loads((HERE/'manifest.json').read_text());M2=json.loads((HERE/'refill_manifest.json').read_text())
 for d in a.dimensions or sorted(map(int,M2)):
  P=Problem(d,a.siblings);assert np.array_equal(P.g,G)
  h=M[str(d)];s=M2[str(d)]
  q0=(HERE/f'global{d}.json.gz').read_bytes();q1=(HERE/f'refill{d}.json.gz').read_bytes()
  assert SHA(q0)==h['witness_sha256'] and SHA(q1)==s['refill_sha256']
  assert SHA(BF)==h['bank_sha256'] and P.ph==h['pointwise_sha256'] and P.oh==h['signed_sha256']
  w=json.loads(gzip.decompress(q0));ff=json.loads(gzip.decompress(q1))
  assert w['dimension']==ff['dimension']==d and w['prior']==PREV[d-50]
  assert w['bank_sha']==SHA(BF) and w['pointwise_sha']==P.ph and w['signed_sha']==P.oh
  if ff['global_sha']!=SHA(q0):
   assert d==60
   prev=(HERE/'global60_cube.json.gz').read_bytes()
   assert SHA(prev)==s['cube_witness_sha256']==ff['global_sha']
   v=json.loads(gzip.decompress(prev))
   assert {k:t for k,t in v.items()if k not in ('lower_bound','tail_gain')}=={k:t for k,t in w.items()if k not in ('lower_bound','tail_gain')}
   assert v['tail_gain']==3583 and w['tail_gain']==4799
   assert ff['final']+1216==w['lower_bound']+ff['gain']
  else:assert ff['final']==w['lower_bound']+ff['gain']
  H,g=advanced() if d==60 else tails(d,a.tail);assert w['tail_gain']==g
  alt,extra=overrides.get(d,(np.empty((0,48),'i8'),0))
  assert extra>=0 and np.all((ip(alt,G)*alt).sum(axis=1)==6)
  assert len(np.unique(key(alt)))==len(alt)
  Hk=np.unique(key(np.r_[H,alt,-alt]));assert np.all(np.sum(ip(H,G)*H,axis=1)==6)
  D={x['class_index']:x['insert'][:]for x in w['changes']}
  assert len(D)==len(w['changes'])
  B0={tuple(x):j for j,x in enumerate(A)};E=[]
  count=0;claimed_extra_retires=0
  for c in ff['changes']:
   i=c['class_index'];assert type(i)is int and 0<=i<len(P.T)
   C=D.setdefault(i,[])
   for L in c['extra']:
    assert len(L)==48 and all(type(z)is int and abs(z)<2**15 for z in L)
    p=tuple(L)
    if p not in B0:B0[p]=n+len(E);E.append(p)
    C.append(B0[p]);count+=1
   claimed_extra_retires+=c['additional_retire']
  assert count==ff['adds'] and claimed_extra_retires==ff['retire'] and ff['gain']==ff['adds']-ff['retire']
  J0=np.r_[A,np.array(E,'i8').reshape(-1,48)]
  DG=ip(J0,G);assert np.all((DG*J0).sum(1)==6)
  both=[]
  for i,I in D.items():
   assert type(i)is int and 0<=i<len(P.T) and len(I)==len(set(I)) and all(type(j)is int and 0<=j<len(J0) for j in I)
   O=fp(DG[I],J0[I].T)
   assert np.all(O.diagonal()==6) and np.all(O[np.triu_indices(len(I),1)]<=2)
   T=P.T[i];assert np.array_equal(ip(ip(T,G),T.T),G)
   both.append(key(ip(J0[I],T)))
  W=np.sort(np.concatenate(both))
  assert np.all(W[1:]!=W[:-1])
  if len(Hk):assert not matched(Hk,W)[1].any()
  assert len(W)==w['adds']+ff['adds']
  RB=[]
  for j in range(0,len(J0),256):RB.extend(np.flatnonzero(row) for row in fp(DG[j:j+256],Ref.T)>2)
  del DG
  old=retire=0;pairs=0;moved_total=0
  for i,T in enumerate(P.T):
   V=P.rows(i);old+=len(V)
   vk=key(ip(V,T));same=matched(W,vk)[1];mark=same.copy()
   moved_total+=int(same.sum())
   if len(Hk):assert not matched(Hk,vk)[1].any()
   if i in D:
    L=D[i];q,isin=matched(rk,key(V));keep=np.flatnonzero(isin)
    remap=np.full(len(Ref),-1,'i8');remap[r0[q[keep]]]=keep
    for j in L:
     Q=remap[RB[j]];mark[Q[Q>=0]]=1
    other=np.flatnonzero(~isin)
    if len(other):mark[other]|=(fp(ip(J0[L],G),V[other].T)>2).any(0)
    pairs+=len(L)*int((~mark).sum())+len(L)*(len(L)-1)//2
   retire+=int(mark.sum())
  assert w['prior']==52416000+POLE[d-50]+old
  assert retire==w['retire']+ff['retire']
  val=w['prior']+len(W)-retire+g
  assert val==w['lower_bound']+ff['gain']
  report=dict(status='PASS',dimension=d,prior=w['prior'],global_gain=w['owner_gain'],refill_gain=ff['gain'],
    tail_gain=g+extra,lower_bound=val+extra,reassigned=moved_total,net_additions=len(W)-retire,added=len(W),retired=retire,exact_same_class_new_pairs=pairs,
    migration_sha=SHA(q0),refill_sha=SHA(q1))
  if expected is not None:
   assert report['lower_bound']==expected[str(d)]
  print(json.dumps(report),flush=True)
if __name__=='__main__':
 if not __debug__:raise SystemExit('assertions required')
 p=argparse.ArgumentParser()
 p.add_argument('--siblings',type=Path,default=HERE.parent)
 p.add_argument('--tail',type=Path,default=HERE.parent/'dim54-60-p48-tail-antipodes')
 p.add_argument('dimensions',type=int,nargs='*')
 run(p.parse_args())
