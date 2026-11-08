"""All 134596 core six-supports, not just supports present in the incumbent."""
from pathlib import Path
from itertools import combinations,product
import numpy as np,argparse,json,time,hashlib,sys
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'common'));from golay import golay24
from optimize_rational import optimize
p=argparse.ArgumentParser();p.add_argument('dimension',type=int);p.add_argument('--cutoff',type=int,default=2);p.add_argument('--seconds',type=int,default=60);p.add_argument('--resume-log',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);st=time.monotonic();d=a.dimension
source=next((ROOT/('verifications/improved/dim27-31-'+('algebraic-completion' if d in (27,31) else 'joint-refills')+'/data')).glob(f'coordinates{d}_*.npz'))
X=np.load(source)['points'].astype(np.int16);cid=np.flatnonzero(~np.any(X[:,24:],axis=1));aid=np.flatnonzero(np.any(X[:,24:],axis=1));core=X[cid,:24];isoc=np.all((core==0)|(abs(core)==4),axis=1)&((core!=0).sum(axis=1)==8)
assert np.max(np.sort(abs(core[~isoc]),axis=1)[:,-6:].sum(axis=1))*4<=64
G=golay24();octads=[tuple(np.flatnonzero(w).tolist()) for w in G[G.sum(axis=1)==8]];assert len(octads)==759
byoc={o:[] for o in octads}
for ix,row in zip(cid[isoc],core[isoc]):byoc[tuple(np.flatnonzero(row).tolist())].append(int(ix))
lookup5={};lookup6={};tab5={};tab6={}
for oi,o in enumerate(octads):
 ids=np.array(byoc[o],np.int64)
 for w in (5,6):
  lookup=lookup5 if w==5 else lookup6;tab=tab5 if w==5 else tab6
  for s in combinations(o,w):
   bits=np.sum((X[ids][:,s]<0).astype(np.uint16)*(1<<np.arange(w)),axis=1,dtype=np.uint16) if len(ids) else np.array([],np.uint16)
   lookup[s]=oi;tab[s]=np.bincount(bits,minlength=1<<w).astype(np.int16)
words=np.arange(64);remove=np.array([(words&((1<<i)-1))|((words>>(i+1))<<i) for i in range(6)])
bodies=[];bb=[];support_count=0;nonseed=0
known_supports={tuple(np.flatnonzero(row[:24]).tolist()) for row in X if np.count_nonzero(row[:24])==6 and np.all((row[:24]==0)|(abs(row[:24])==4))}
for s in combinations(range(24),6):
 if s in tab6:cnt=tab6[s];related=[lookup6[s]]
 else:
  five=[s[:i]+s[i+1:] for i in range(6)];cnt=sum(tab5[q][remove[i]] for i,q in enumerate(five));related=[lookup5[q] for q in five]
 good=np.flatnonzero(cnt<=a.cutoff)
 if not len(good):continue
 support_count+=1;nonseed+=int(s not in known_supports);ids=np.array([j for o in set(related) for j in byoc[octads[o]]],np.int64)
 for word in good:
  row=np.zeros(d,np.int16);row[list(s)]=np.where((word>>np.arange(6))&1,-4,4)
  # Exact local oracle equality checks all potentially violating core octads.
  blocked=ids[X[ids].astype(np.int64)@row.astype(np.int64)>64].tolist();assert len(blocked)==cnt[word]
  bodies.append(row);bb.append(blocked)
print(json.dumps({'all_supports':134596,'admitted_supports':support_count,'admitted_nonseed_supports':nonseed,'admitted_signed_bodies':len(bodies),'seconds':time.monotonic()-st}),flush=True)
known={row.tobytes() for row in X};tail=[]
for i,j in combinations(range(24,d),2):
 for s,t in product((-4,4),repeat=2):
  v=np.zeros(d,np.int16);v[i]=s;v[j]=t;tail.append(v)
B=np.array(bodies,np.int16);T=np.array(tail,np.int16);rng=np.random.default_rng(20261108);order=rng.permutation(aid);pool=[];bs=[];records=[];checks=0
skip=set()
if a.resume_log:
 for line in a.resume_log.read_text().splitlines():
  try: rec=json.loads(line)
  except json.JSONDecodeError: continue
  if 'body_start' in rec and rec['new_admitted']==0:skip.add(rec['body_start'])
for at in range(0,len(B),4096):
 if at in skip:continue
 Y=(B[at:at+4096,None,:]+T[None,:,:]).reshape((-1,d));nc=min(4096,len(B)-at);counts=np.repeat(np.array([len(b) for b in bb[at:at+4096]],np.int16)[:,None],len(T),axis=1);active_bodies=np.arange(nc)
 for off in range(0,len(order),256):
  if not len(active_bodies):break
  oldrows=X[order[off:off+256]];bp=B[at+active_bodies,:24].astype(np.float32)@oldrows[:,:24].astype(np.float32).T;td=T[:,24:].astype(np.float32)@oldrows[:,24:].astype(np.float32).T
  inc=np.zeros((len(active_bodies),len(T)),np.float32)
  for value in np.unique(td):
   inc+=(bp>64-value).astype(np.float32)@(td==value).astype(np.float32).T
  assert np.all(inc==np.rint(inc))
  counts[active_bodies]+=inc.astype(np.int16);checks+=len(active_bodies)*len(T)*len(oldrows);counts=np.minimum(counts,a.cutoff+1)
  active_bodies=active_bodies[np.any(counts[active_bodies]<=a.cutoff,axis=1)]
 active=np.flatnonzero(counts.ravel()<=a.cutoff)
 new=[i for i in active if Y[i].tobytes() not in known]
 for row in Y[new]:
  dot=np.zeros(len(X),np.int64)
  for j in np.flatnonzero(row):dot+=int(row[j])*X[:,j].astype(np.int64)
  b=np.flatnonzero(dot>64).tolist();assert len(b)<=a.cutoff;pool.append(row);bs.append(b)
 rec={'body_start':at,'tested':len(Y),'new_admitted':len(new),'seconds':time.monotonic()-st};records.append(rec);print(json.dumps(rec),flush=True)
 np.savez_compressed(a.out/'working_pool.npz',candidates=np.array(pool,np.int16).reshape((-1,d)),blocker_indptr=np.r_[0,np.cumsum(list(map(len,bs)))],blocker_indices=np.array([j for b in bs for j in b],np.int64))
C=np.array(pool,np.int64).reshape((-1,d));np.savez_compressed(a.out/'pool.npz',candidates=C,blocker_indptr=np.r_[0,np.cumsum(list(map(len,bs)))],blocker_indices=np.array([j for b in bs for j in b],np.int64));print('OPTIMIZE',len(C),flush=True)
r=optimize(C,bs,128,a.seconds) if len(C) else {'gain':0}
if r.get('gain',0)>0:
 keep=np.ones(len(X),bool);keep[r['removed']]=False;out=np.r_[X[keep],C[r['added']]].astype(np.int16);np.savez_compressed(a.out/'candidate.npz',points=out,squared_norm=128,inner_product_limit=64,coordinate_denominator=2)
rep={'dimension':d,'source':str(source.relative_to(ROOT)),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'baseline':len(X),'all_six_supports':134596,'admitted_supports':support_count,'nonseed_supports':nonseed,'signed_bodies':len(B),'tails':len(T),'raw_candidates':len(B)*len(T),'admitted':len(C),'reused_zero_chunks':len(skip),'resume_log':str(a.resume_log) if a.resume_log else None,'cutoff':a.cutoff,'checks':checks,'optimization':r,'seconds':time.monotonic()-st,'scope':'Every core six-support and every sign, with stated core/global blocker cutoff; all signed tail roots. No support selected from source required. Candidate blockers independently rebuilt by full int64 sums.'};(a.out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep,indent=2),flush=True)
