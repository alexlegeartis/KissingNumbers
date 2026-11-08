"""Complete negative-contact closure anchored on an entire latitude layer."""
from pathlib import Path
import numpy as np,argparse,json,time,hashlib,importlib.util
ROOT=Path(__file__).resolve().parents[2]
p=argparse.ArgumentParser();p.add_argument('dimension',type=int);p.add_argument('--tail-norm',type=int,default=64);p.add_argument('--source',type=Path);p.add_argument('--unpaired',action='store_true');p.add_argument('--cutoff',type=int,default=2);p.add_argument('--out',type=Path,required=True);p.add_argument('--seconds',type=float,default=90);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);st=time.monotonic();d=a.dimension
source=a.source or next((ROOT/('verifications/improved/dim27-31-'+('algebraic-completion' if d in (27,31) else 'joint-refills')+'/data')).glob(f'coordinates{d}_*.npz'))
X=np.load(source)['points'].astype(np.int16);F=X.astype(np.float32);T=np.sum(X[:,24:].astype(np.int64)**2,axis=1);known={v.tobytes() for v in X};anchors=np.array([i for i,row in enumerate(X) if (-row).tobytes() not in known]) if a.unpaired else np.flatnonzero(T==a.tail_norm);seen=set();pool=[];blockers=[];rng=np.random.default_rng(20261108);order=rng.permutation(len(X));records=[];generated=0;checks=0
# Every X_i+X_j with product -64 has norm128 exactly. All signs included.
for at in range(0,len(anchors),32):
 ids=anchors[at:at+32];R=F[ids]@F.T;ii,jj=np.nonzero(R==-64);Y=X[ids[ii]]+X[jj];generated+=len(Y)
 new=[]
 for row in Y:
  key=row.tobytes()
  if key not in known and key not in seen:seen.add(key);new.append(row)
 Y=np.array(new,np.int16).reshape((-1,d));active=np.arange(len(Y));counts=np.zeros(len(Y),np.int16)
 assert np.all(np.sum(Y.astype(np.int64)**2,axis=1)==128)
 for off in range(0,len(order),256):
  if not len(active):break
  bad=Y[active].astype(np.float32)@F[order[off:off+256]].T>64;checks+=bad.size;counts[active]+=bad.sum(axis=1).astype(np.int16);active=active[counts[active]<=a.cutoff]
 bs=[]
 for row in Y[active]:
  dot=np.zeros(len(X),np.int64)
  for j in np.flatnonzero(row):dot+=int(row[j])*X[:,j].astype(np.int64)
  b=np.flatnonzero(dot>64).tolist();assert len(b)<=a.cutoff;bs.append(b)
 pool.extend(Y[active]);blockers.extend(bs);rec={'anchor_start':at,'anchors':len(ids),'sum_incidences':len(ii),'new_sums':len(Y),'admitted':len(active),'free':sum(not b for b in bs),'seconds':time.monotonic()-st};records.append(rec);print(json.dumps(rec),flush=True)
 (a.out/'partial.json').write_text(json.dumps(records,indent=2)+'\n')
C=np.array(pool,np.int64).reshape((-1,d));np.savez_compressed(a.out/'pool.npz',candidates=C,blocker_indptr=np.r_[0,np.cumsum(list(map(len,blockers)))],blocker_indices=np.array([j for b in blockers for j in b],np.int64));print('OPTIMIZE',len(C),flush=True)
spec=importlib.util.spec_from_file_location('cc',ROOT/'research/sol_discovery_20261104/contact_campaign.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
from optimize_rational import optimize
r=optimize(C,blockers,128,a.seconds) if len(C) else {'gain':0}
if r.get('gain',0)>0:
 keep=np.ones(len(X),bool);keep[r['removed']]=False;out=np.r_[X[keep],C[r['added']]].astype(np.int16);np.savez_compressed(a.out/'candidate.npz',points=out,squared_norm=128,inner_product_limit=64,coordinate_denominator=2)
rep={'dimension':d,'source':str(source.relative_to(ROOT)),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'baseline':len(X),'tail_norm':None if a.unpaired else a.tail_norm,'unpaired_anchor_family':a.unpaired,'anchor_count':len(anchors),'generated_incidences':generated,'unique_new_sums':len(seen),'admitted':len(C),'cutoff':a.cutoff,'dot_products':checks,'optimization':r,'seconds':time.monotonic()-st,'scope':'Complete sum closure for every incumbent point in stated latitude layer paired with all incumbent points; blocker cutoff restricts optimization pool; numerical non-gain is not a global bound.'};(a.out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep,indent=2),flush=True)
