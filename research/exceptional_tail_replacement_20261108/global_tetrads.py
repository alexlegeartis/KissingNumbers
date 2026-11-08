"""Whole norm64-body/norm64-tail signed tetrad grammar, all 24 core axes."""
from pathlib import Path
from itertools import combinations,product
import numpy as np,argparse,json,hashlib,time,importlib.util
ROOT=Path(__file__).resolve().parents[2]
p=argparse.ArgumentParser();p.add_argument('dimension',type=int);p.add_argument('--cutoff',type=int,default=4);p.add_argument('--seconds',type=int,default=120);p.add_argument('--start-tail',type=int,default=0);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);st=time.monotonic();d=a.dimension
source=next((ROOT/('verifications/improved/dim27-31-'+('algebraic-completion' if d in (27,31) else 'joint-refills')+'/data')).glob(f'coordinates{d}_*.npz'))
X=np.load(source)['points'].astype(np.int16);core=X[~np.any(X[:,24:],axis=1),:24];aid=np.flatnonzero(np.any(X[:,24:],axis=1));known={v.tobytes() for v in X};rng=np.random.default_rng(20261108);order=rng.permutation(aid)
# For every core point, no signed tetrad can have product above 64.
assert np.max(np.sort(abs(core),axis=1)[:,-4:].sum(axis=1))*4<=64
signs=np.array(list(product((-4,4),repeat=4)),np.int16)
B=[]
for support in combinations(range(24),4):
 z=np.zeros((16,d),np.int16);z[:,support]=signs;B.append(z)
B=np.concatenate(B);tails=[]
for support in combinations(range(24,d),4):
 z=np.zeros((16,d),np.int16);z[:,support]=signs;tails.extend(z)
# Include the other universally compatible norm64 tail shape: one coordinate 8.
for i in range(24,d):
 for s in (-8,8):
  z=np.zeros(d,np.int16);z[i]=s;tails.append(z)
raw=len(B)*len(tails);Ys=[];blockers=[];records=[];checks=0
for ti in range(a.start_tail,len(tails),8):
 T=np.array(tails[ti:ti+8]);Y=(B[:,None,:]+T[None,:,:]).reshape((-1,d));ids=np.arange(len(Y));counts=np.zeros(len(Y),np.int16)
 for at in range(0,len(order),256):
  if not len(ids):break
  active=Y[ids].astype(np.float32);bad=active@X[order[at:at+256]].astype(np.float32).T>64;checks+=bad.size
  counts[ids]+=bad.sum(axis=1).astype(np.int16);ids=ids[counts[ids]<=a.cutoff]
 new=[j for j in ids if Y[j].tobytes() not in known];Z=Y[new]
 # Independently reconstruct complete blocker lists for the admitted pool.
 bs=[]
 for row in Z:
  dot=np.zeros(len(X),np.int64)
  for j in np.flatnonzero(row):dot+=int(row[j])*X[:,j].astype(np.int64)
  b=np.flatnonzero(dot>64).tolist();assert len(b)<=a.cutoff;bs.append(b)
 Ys.extend(Z);blockers.extend(bs);rec={'tail_start':ti,'tails':len(T),'tested':len(Y),'retained':len(ids),'new':len(Z),'zero_blockers':sum(not b for b in bs),'seconds':time.monotonic()-st};records.append(rec);print(json.dumps(rec),flush=True)
 (a.out/'partial.json').write_text(json.dumps(records,indent=2)+'\n')
 np.savez_compressed(a.out/'working_pool.npz',candidates=np.array(Ys,np.int16).reshape((-1,d)),blocker_indptr=np.r_[0,np.cumsum(list(map(len,blockers)))],blocker_indices=np.array([j for b in blockers for j in b],np.int64))
Z=np.array(Ys,np.int64).reshape((-1,d));np.savez_compressed(a.out/'pool.npz',candidates=Z,blocker_indptr=np.r_[0,np.cumsum(list(map(len,blockers)))],blocker_indices=np.array([j for b in blockers for j in b],np.int64))
print('OPTIMIZE',len(Z),flush=True)
spec=importlib.util.spec_from_file_location('contact_opt',ROOT/'research/sol_discovery_20261104/contact_campaign.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
from optimize_rational import optimize
r=optimize(Z,blockers,128,a.seconds) if len(Z) else {'gain':0}
if r.get('gain',0)>0:
 keep=np.ones(len(X),bool);keep[r['removed']]=False;F=np.r_[X[keep],Z[r['added']]].astype(np.int16);np.savez_compressed(a.out/'candidate.npz',points=F,squared_norm=128,inner_product_limit=64,coordinate_denominator=2)
rep={'start_tail':a.start_tail,'dimension':d,'source':str(source.relative_to(ROOT)),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'baseline':len(X),'core_tetrads':10626,'tail_vectors':len(tails),'raw_candidates':raw,'new_admitted':len(Z),'zero_blockers':sum(not b for b in blockers),'cutoff':a.cutoff,'products':checks,'optimization':r,'seconds':time.monotonic()-st,'scope':'All signed four-core/four-tail and four-core/one-heavy-tail norm128 points. All incumbent points deletable; blocker cutoff bounds the finite pool only; no unrestricted exclusion.'};(a.out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep,indent=2),flush=True)
