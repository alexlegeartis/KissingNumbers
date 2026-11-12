"""Bounded individual-point exchange search, with exact finite conflict data.

Pool selection and annealing are heuristics. No negative outcome is a proof.
Default core triples comprise the mixed {11,23,j} family and twenty additional
mixed triples selected by minimum complete-sign-bundle deletion cost.
"""
import argparse,itertools,json,time,subprocess,hashlib
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from source_paths import PACKINGSTAR_DEFAULT,packingstar_file

EXTRA=[(2,9,10),(2,10,15),(2,7,10),(2,10,13),(2,3,10),(2,10,19),(2,10,17),(0,17,23),(0,3,23),(2,10,16),(0,9,23),(0,7,23),(2,10,20),(2,8,10),(2,10,14),(2,10,12),(2,10,18),(2,10,21),(2,10,22),(0,21,23)]
ap=argparse.ArgumentParser();ap.add_argument('--packingstar',type=Path,default=PACKINGSTAR_DEFAULT);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--report',type=Path,required=True);ap.add_argument('--wall-seconds',type=float,default=290);ap.add_argument('--seed',type=int,default=20261112);ap.add_argument('--search-mode',choices=['anneal','continuation'],default='continuation');args=ap.parse_args()
start=time.monotonic();args.output.mkdir(parents=True,exist_ok=True);args.report.parent.mkdir(parents=True,exist_ok=True)
source=packingstar_file(args.packingstar,31);f=np.load(source,allow_pickle=False);X=f['points'].astype(np.int16);assert int(f['coordinate_denominator'])==2 and np.all((X.astype(np.int32)**2).sum(1)==128)
report={'status':'running','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'source_points':len(X),'seed':args.seed,'search_mode':args.search_mode,'wall_budget_seconds':args.wall_seconds,'proof_scope':'heuristic search; negative outcomes are unresolved, not no-gain certificates','core_triples':sorted(set(tuple(sorted((11,23,j))) for j in range(24) if j not in [11,23])|set(tuple(sorted(t)) for t in EXTRA))}
def checkpoint(stage,**kw):
 report.update(stage=stage,elapsed_seconds=time.monotonic()-start,**kw);args.report.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'stage':stage,**kw,'elapsed':report['elapsed_seconds']}),flush=True)
def famdata(r):
 supports=report['core_triples'] if r==3 else list(itertools.combinations(range(24),r))
 H=[]
 for hs in supports:
  for ss in itertools.product([-1,1],repeat=r):
   h=np.zeros(24,np.int8);h[list(hs)]=ss;H.append(h)
 T=[]
 for ts in itertools.combinations(range(7),8-r):
  for ss in itertools.product([-1,1],repeat=8-r):
   t=np.zeros(7,np.int8);t[list(ts)]=ss;T.append(t)
 H=np.array(H);T=np.array(T);am=np.sort(np.abs(X[:,:24]),axis=1)[:,-r:].sum(1)+np.sort(np.abs(X[:,24:]),axis=1)[:,-(8-r):].sum(1)
 old=np.flatnonzero(am>16);tailmax=np.sort(np.abs(X[old,24:]),axis=1)[:,-(8-r):].sum(1)
 td=T.astype(np.float32)@X[old,24:].astype(np.float32).T;hd=H.astype(np.float32)@X[old,:24].astype(np.float32).T
 costs=np.empty((len(H),len(T)),np.int32)
 for i in range(len(H)):
  a=np.flatnonzero(hd[i]+tailmax>16);costs[i]=(td[:,a]>16-hd[i,a]).sum(1)
 return {'r':r,'H':H,'T':T,'old':old,'tailmax':tailmax,'td':td,'hd':hd,'costs':costs}
families=[]
for r in [3,2,1]:
 F=famdata(r);families.append(F);checkpoint('priced_family',family=r,priced_points=sum(z['costs'].size for z in families),family_points=F['costs'].size,family_minimum=int(F['costs'].min()),family_median=float(np.median(F['costs'])),family_maximum=int(F['costs'].max()))
report['family_costs']=[{'r':z['r'],'points':z['costs'].size,'minimum':int(z['costs'].min()),'median':float(np.median(z['costs'])),'maximum':int(z['costs'].max())} for z in families]
def blockers(F,h,t):
 a=np.flatnonzero(F['hd'][h]+F['tailmax']>16)
 return F['old'][a[F['td'][t,a]>16-F['hd'][h,a]]]
# Retain alternatives for every signed tail; then adapt using old-blocker reuse.
seed=[]
for fi,F in enumerate(families):
 order=np.argsort(F['costs'],axis=0,kind='stable')[:8]
 for t in range(len(F['T'])):
  for h in order[:,t]:seed.append((fi,int(h),t))
freq=np.zeros(len(X),np.int32)
for fi,h,t in seed:freq[blockers(families[fi],h,t)]+=1
weights=1/np.sqrt(1+freq.astype(float));pool=set(seed)
for fi,F in enumerate(families):
 scores=np.empty(F['costs'].shape,np.float32)
 for h in range(len(F['H'])):
  a=np.flatnonzero(F['hd'][h]+F['tailmax']>16)
  bad=F['td'][:,a]>16-F['hd'][h,a]
  scores[h]=bad.astype(np.float32)@weights[F['old'][a]].astype(np.float32)
 extra={3:8000,2:4000,1:2000}[F['r']]
 for z in np.argsort(scores.ravel(),kind='stable')[:extra]:
  h,t=divmod(int(z),len(F['T']));pool.add((fi,h,t))
pool=sorted(pool);Y=[];B=[]
for fi,h,t in pool:
 F=families[fi];Y.append(np.r_[F['H'][h],F['T'][t]]*4);B.append(blockers(F,h,t).astype(np.int32))
Y=np.array(Y,np.int8);bp=np.r_[0,np.cumsum([len(b) for b in B])].astype(np.int64);bi=np.concatenate(B);victims,bi=np.unique(bi,return_inverse=True);bi=bi.astype(np.int32)
checkpoint('adapted_pool',pool_points=len(Y),pool_old_victims=len(victims),pool_old_incidences=len(bi),seed_points=len(seed))
# Binary32 products are exact here: integral operands and every partial sum
# has absolute value at most 128, well below the 2**24 exact-integer range.
rows=[];cols=[];yf=Y.astype(np.float32)
for at in range(0,len(Y),256):
 bad=yf[at:at+256]@yf.T>64
 ii,jj=np.nonzero(bad);ok=at+ii!=jj;rows.append((at+ii[ok]).astype(np.int32));cols.append(jj[ok].astype(np.int32))
rr=np.concatenate(rows);cc=np.concatenate(cols);A=coo_matrix((np.ones(len(rr),np.uint8),(rr,cc)),shape=(len(Y),len(Y))).tocsr()
checkpoint('conflict_graph',new_conflict_edges=A.nnz//2)
binary=args.output/'model.bin'
with binary.open('wb') as out:
 np.array([len(Y),len(victims),len(bi),A.nnz],np.int64).tofile(out);bp.tofile(out);bi.tofile(out);A.indptr.astype(np.int64).tofile(out);A.indices.astype(np.int32).tofile(out)
np.savez_compressed(args.output/'pool.npz',points=Y,old_indices=victims,blocker_indptr=bp,blocker_indices=bi,metadata=np.array(pool,np.int32))
worker=args.output/'partial_sign_search';subprocess.run(['g++','-O3','-std=c++17',str(Path(__file__).with_name('partial_sign_search.cpp')),'-o',str(worker)],check=True,timeout=20)
remaining=max(1,args.wall_seconds-(time.monotonic()-start)-8);result=args.output/'search_result.json'
try:
 subprocess.run([str(worker),str(binary),str(remaining),str(args.seed),str(result),args.search_mode],check=True,timeout=remaining+3)
except subprocess.TimeoutExpired:report['external_search_timeout']=True
if result.exists():
 s=json.loads(result.read_text());chosen=np.array(s.pop('selected'),np.int32);near=np.array(s.pop('best_nonempty_selected'),np.int32);report['search']=s
 for label,selection in [('best',chosen),('best_nonempty',near)]:
  removed=np.unique(np.concatenate([B[i] for i in selection])) if len(selection) else np.empty(0,np.int32)
  assert len(selection)-len(removed)==(s['gain'] if label=='best' else s['best_nonempty_gain'])
  z=Y[selection].astype(np.int32);assert len(z)==0 or np.all((z@z.T-np.eye(len(z),dtype=np.int32)*128)<=64)
  report[label]={'additions':len(selection),'deletions':len(removed),'gain':len(selection)-len(removed),'new_new_exact_replay':True,'all_blockers_removable':True}
  if label=='best' and len(selection)-len(removed)>0:
   keep=np.ones(len(X),bool);keep[removed]=False;new=np.r_[X[keep],Y[selection]].astype(np.int16)
   for at in range(0,len(selection),256):
    dots=Y[selection[at:at+256]].astype(np.float32)@new.astype(np.float32).T
    for j in range(len(dots)):dots[j,len(X)-len(removed)+at+j]=0
    assert np.max(dots)<=64
   np.savez_compressed(args.output/f'candidate31_{len(new)}.npz',points=new,squared_norm=128,inner_product_limit=64,coordinate_denominator=2)
   report['positive_changed_pair_replay']='PASS; integer-exact binary32 products with bounded integer partial sums'
 checkpoint('finished',status='heuristic completed; no global conclusion',best_gain=s['gain'])
else:checkpoint('finished',status='unresolved: no heuristic checkpoint returned')
