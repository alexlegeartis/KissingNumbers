"""Complete one-missing-parent Leech midpoint catalogue, with an orbit cost cap."""
import argparse,itertools,time,json
from pathlib import Path
import numpy as np
from image_search import above

def main():
 ap=argparse.ArgumentParser();ap.add_argument('dim',type=int);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--budget',type=int,default=2);ap.add_argument('--limit',type=int,default=1000000);a=ap.parse_args();a.out.mkdir(exist_ok=False);r=(Path(__file__).resolve().parent/'data');start=time.monotonic();d=a.dim;baseline=201180 if d==26 else 209862;z=np.load(r/f'PackingStar/25D-31D/new_bounds_20261111/D{d}_{baseline}_coordinates.npz')
 if d==26:O=2*np.r_[z['rational_points'],z['quadratic_A']];OB=2*np.r_[np.zeros_like(z['rational_points']),z['quadratic_B']]
 else:O=2*z['points'];OB=np.zeros_like(O)
 L=np.load(r/'leech.npy');g=np.load(r/f'parentgraph{d}/graph.npz');sf=g['source_for_leech'];lk={v.tobytes():i for i,v in enumerate(L)};positive=L[np.arange(len(L)),np.argmax(L!=0,axis=1)]>0;canon=np.array([i if positive[i] else lk[(-v).tobytes()] for i,v in enumerate(L)]);ret=np.flatnonzero(sf>=0);us=np.flatnonzero((sf<0)&positive);deg=np.bincount(g['parents'].ravel(),minlength=len(L));us=us[np.argsort(deg[us],kind='stable')];LF=L[ret].astype(np.float32);oldH=g['bodies'];pos=oldH[np.arange(len(oldH)),np.argmax(oldH!=0,axis=1)]>0;oldH=oldH[pos];oldT=g['tails'][pos];T=[]
 for i,j in itertools.combinations(range(d-24),2):
  for s in [-1,1]:
   v=np.zeros(d-24,np.int16);v[[i,j]]=[4,4*s];T.append(v)
 T=np.array(T);tl={v.tobytes():i for i,v in enumerate(T)};oldcolor=[]
 for t in oldT:
  t=t.copy()
  if t[np.flatnonzero(t)[0]]<0:t=-t
  oldcolor.append(tl[t.tobytes()])
 oldcolor=np.array(oldcolor);cost=np.abs(T.astype(np.int64)@T.T);core=np.sum(O[:,:24].astype(np.int64)**2+3*OB[:,:24].astype(np.int64)**2,axis=1);hg={};qg={}
 for i in np.flatnonzero(core==64):
  h=O[i,:24].copy();t=O[i,24:].copy()
  if h[np.flatnonzero(h)[0]]<0:h=-h
  if t[np.flatnonzero(t)[0]]<0:t=-t
  hg[np.r_[h,t].tobytes()]=(h,t)
 for i in np.flatnonzero(np.any(OB,axis=1)):
  h=OB[i,:24].copy();t=O[i,24:].copy()
  if h[np.flatnonzero(h)[0]]<0:h=-h
  if t[np.flatnonzero(t)[0]]<0:t=-t
  qg[np.r_[h,t].tobytes()]=(h,t)
 HA=np.array([v[0] for v in hg.values()],np.float32);HT=np.array([v[1] for v in hg.values()],np.float32);QB=np.array([v[0] for v in qg.values()],np.float32).reshape(-1,24);QT=np.array([v[1] for v in qg.values()],np.float32).reshape(-1,d-24);hc=np.abs(T.astype(np.float32)@HT.T);qc=np.abs(T.astype(np.float32)@QT.T);order=np.random.default_rng(349).permutation(len(oldH));oldH=oldH[order].astype(np.float32);oldcolor=oldcolor[order];poolA=[];poolP=[];poolcost=[];tested=0
 for number,u in enumerate(us[:a.limit]):
  ps=ret[LF@L[u].astype(np.float32)==16];H=L[ps]+L[u];tested+=len(H);head=np.abs(H.astype(np.float32)@HA.T);counts=np.array([(head+hc[c]>64).sum(1) for c in range(len(T))],np.int16).T
  if len(QB):
   qu=np.abs(H.astype(np.float32)@QB.T)
   for c in range(len(T)):counts[:,c]+=above(qc[c][None,:],qu).sum(1)
  active=np.flatnonzero((counts<=a.budget).any(1))
  for at in range(0,len(oldH),128):
   if not len(active):break
   D=np.abs(H[active].astype(np.float32)@oldH[at:at+128].T);co=oldcolor[at:at+128]
   for c in range(len(T)):counts[active,c]+=(D+cost[c,co][None,:]>64).sum(1)
   active=active[(counts[active]<=a.budget).any(1)]
  for j,c in np.argwhere(counts<=a.budget):
   if j not in active:continue
   h=H[j].copy()
   if h[np.flatnonzero(h)[0]]<0:h=-h
   ca=np.r_[h,T[c]].astype(np.int16);poolA.append(ca);poolP.append(np.sort(canon[[u,ps[j]]]));poolcost.append(int(counts[j,c]))
   if counts[j,c]==0:
    remove=np.array([int(sf[ps[j]]),int(sf[lk[(-L[ps[j]]).tobytes()]])]);keep=np.ones(len(O),bool);keep[remove]=False;C=np.array([np.r_[s*h,t*T[c]] for s,t in itertools.product([-1,1],repeat=2)],np.int16);RA=np.r_[O[keep],C];RB=np.r_[OB[keep],np.zeros_like(C)]
    for z,v in enumerate(C):
     fail=above(RA.astype(np.int64)@v.astype(np.int64),RB.astype(np.int64)@v.astype(np.int64));fail[len(O)-len(remove)+z]=False;assert not fail.any()
    np.savez_compressed(a.out/'candidate.npz',rational=RA,radical=RB,radicand=3,squared_norm=128,inner_product_limit=64,coordinate_denominator=2);np.savez_compressed(a.out/'exchange.npz',added_A=C,added_B=np.zeros_like(C),removed_indices=remove);print('POSITIVE_EXACT_CHANGED_PAIRS',len(RA),flush=True);return
  if number%25==0:print('PROGRESS',number,len(us),'bodies',tested,'states',len(poolA),'seconds',round(time.monotonic()-start,1),flush=True)
 np.savez_compressed(a.out/'pool.npz',candidate_A=np.array(poolA,np.int16).reshape(-1,d),candidate_B=np.zeros((len(poolA),d),np.int16),parent_indices=np.array(poolP,np.int32).reshape(-1,2),orbit_blocker_counts=poolcost,radicand=1);report=dict(dimension=d,missing_parents=min(len(us),a.limit),bodies_tested=tested,states=len(poolA),budget=a.budget,histogram={int(k):int(v) for k,v in zip(*np.unique(poolcost,return_counts=True))},seconds=time.monotonic()-start);(a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('FINISHED',json.dumps(report),flush=True)
if __name__=='__main__':main()
