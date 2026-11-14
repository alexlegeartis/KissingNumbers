"""Price complete contact cycles against the incumbent and seek one-for-two trades."""
import argparse,itertools,time,json
from pathlib import Path
from collections import defaultdict
import numpy as np
from image_search import above

def main():
 ap=argparse.ArgumentParser();ap.add_argument('pack',type=Path);ap.add_argument('catalog',type=Path);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--max-blockers',type=int,default=4);a=ap.parse_args();a.out.mkdir(exist_ok=False);start=time.monotonic();r=(Path(__file__).resolve().parent/'data');z=np.load(a.pack/'data.npz');d=int(z['dimension']);T=z['tail_palette'];nt=len(T);L=z['leech'];F=z['fixed_A'];FB=z['fixed_B'];blocked=z['blocked'];flat=z['state'][z['initial']];ci=flat//nt;phase=flat%nt;OH=z['H'].reshape(-1,4,24)[ci];OO=z['offset'].reshape(-1,4)[ci];OP=z['parents'][ci];oldcolors=(OO^phase[:,None]).ravel();oldHF=OH.reshape(-1,24).astype(np.float32);nold=len(ci);lk={v.tobytes():i for i,v in enumerate(L)};anti=np.array([lk[(-v).tobytes()] for v in L]);owner=np.full(len(L),-1,np.int32);owner[OP.ravel()]=np.repeat(np.arange(nold),4);cost=np.abs(T.astype(np.int64)@T.T);c=np.load(a.catalog);HH=c['H'];off=c['offset'];PV=c['parent_vectors'];n=len(HH);parent=np.array([[lk[(v if v[np.flatnonzero(v)[0]]>0 else -v).tobytes()] for v in ps] for ps in PV],np.int32);assert not blocked[parent].any();buckets=defaultdict(lambda:defaultdict(list));kept=[];hist={};pairs_tried=0
 def compatible(x,y):
  if len(set(map(int,x[2]))&set(map(int,y[2]))):return False
  return np.all(np.abs(x[0].astype(np.int64)@y[0].T)+cost[(x[1]^x[3])[:,None],(y[1]^y[3])[None,:]]<=64)
 def positive(records,deleted):
  keepold=~np.isin(np.arange(nold),list(deleted));allH=np.r_[OH[keepold],np.array([v[0] for v in records])];allO=np.r_[OO[keepold],np.array([v[1] for v in records])];allP=np.r_[OP[keepold],np.array([v[2] for v in records])];allPh=np.r_[phase[keepold],np.array([v[3] for v in records])];assert len(np.unique(allP))==4*len(allP);keep=~blocked.copy();keep[allP.ravel()]=False;keep[anti[allP.ravel()]]=False;C=np.array([np.r_[s*h,t*T[int(o)^int(ph)]] for hs,os,ph in zip(allH,allO,allPh) for h,o in zip(hs,os) for s,t in itertools.product([-1,1],repeat=2)],np.int16);LC=np.c_[2*L[keep],np.zeros((int(keep.sum()),d-24),np.int16)];RA=np.r_[LC,F,C];RB=np.r_[np.zeros_like(LC),FB,np.zeros_like(C)];expected=(201180 if d==26 else 209862)+8*(len(records)-len(deleted));assert len(RA)==expected;assert np.all(np.sum(RA.astype(np.int64)**2+3*RB.astype(np.int64)**2,axis=1)==128) and not np.any(np.sum(RA.astype(np.int64)*RB.astype(np.int64),axis=1));src=np.load(r/f'PackingStar/25D-31D/new_bounds_20261111/D{d}_{201180 if d==26 else 209862}_coordinates.npz')
  if d==26:SA=2*np.r_[src['rational_points'],src['quadratic_A']];SB=2*np.r_[np.zeros_like(src['rational_points']),src['quadratic_B']]
  else:SA=2*src['points'];SB=np.zeros_like(SA)
  sourcekeys={np.r_[v,b].tobytes() for v,b in zip(SA,SB)};keys=[np.r_[v,b].tobytes() for v,b in zip(RA,RB)];assert len(set(keys))==len(RA);added=np.flatnonzero([k not in sourcekeys for k in keys]);newkeys=set(keys);removed=np.flatnonzero([np.r_[v,b].tobytes() not in newkeys for v,b in zip(SA,SB)])
  for j in added:
   ca=RA[j].astype(np.int64);cb=RB[j].astype(np.int64);fail=above(RA.astype(np.int64)@ca+3*RB.astype(np.int64)@cb,RA.astype(np.int64)@cb+RB.astype(np.int64)@ca);fail[j]=False;assert not fail.any(),('changed geometry',j,np.flatnonzero(fail)[:8])
  np.savez_compressed(a.out/'candidate.npz',rational=RA,radical=RB,radicand=3,squared_norm=128,inner_product_limit=64,coordinate_denominator=2);np.savez_compressed(a.out/'exchange.npz',added_A=RA[added],added_B=RB[added],removed_indices=removed);report=dict(count=len(RA),gain=len(RA)-len(SA),added=len(added),removed=len(removed),new_cycles=len(records),deleted_cycles=list(map(int,deleted)),seconds=time.monotonic()-start);(a.out/'positive.json').write_text(json.dumps(report,indent=2)+'\n');print('POSITIVE_EXACT_CHANGED_PAIRS',json.dumps(report),flush=True)
 print('LOADED',n,'cycles','baseline',nold,flush=True)
 for at in range(0,n,128):
  H=HH[at:at+128];o=off[at:at+128];P=parent[at:at+128];nb=len(H);D=np.abs(H.reshape(-1,24).astype(np.float32)@oldHF.T).reshape(nb,4,nold,4);pc=np.abs(H.reshape(-1,24).astype(np.float32)@F[:,:24].astype(np.float32).T).reshape(nb,4,-1);qc=np.abs(H.reshape(-1,24).astype(np.float32)@FB[:,:24].astype(np.float32).T).reshape(nb,4,-1);owners=owner[P]
  parity=o[:,:,None,None]^(oldcolors%2).reshape(1,1,nold,4);weak=D>32;weak0=(weak&(parity==0)).any((1,3));weak1=(weak&(parity==1)).any((1,3));strong=(D>48).any((1,3));hard=(D>64).any((1,3));rr,cc=np.nonzero(owners>=0);hard[rr,owners[rr,cc]]=True;oldplane=oldcolors.reshape(nold,4)[:,0]//2
  for ph in range(nt):
   newcol=o^ph;tc=np.abs(T[newcol].astype(np.float32)@F[:,24:].astype(np.float32).T);allow=~above(pc+tc,qc).any((1,2));bad=hard|np.where(oldplane[None,:]==ph//2,weak1 if ph&1 else weak0,strong);counts=bad.sum(1)
   for j in np.flatnonzero(allow&(counts<=a.max_blockers)):
    bs=np.flatnonzero(bad[j]);rec=(H[j].copy(),o[j].copy(),P[j].copy(),int(ph),int(at+j));hist[len(bs)]=hist.get(len(bs),0)+1;kept.append((rec,bs))
    if len(bs)==0:positive([rec],[]);return
    if len(bs)!=1:continue
    b=int(bs[0]);mask=sum(1<<k for k,p in enumerate(OP[b]) if p in P[j])
    for oldmask,previous in buckets[b].items():
     if mask&oldmask:continue
     for prev in previous:
      pairs_tried+=1
      if compatible(rec,prev):positive([prev,rec],[b]);return
    buckets[b][mask].append(rec)
  if at%4096==0:print('PROGRESS',at,n,'kept',len(kept),'hist',hist,'pairs',pairs_tried,'seconds',round(time.monotonic()-start,1),flush=True)
 recs=[v[0] for v in kept];bs=[v[1] for v in kept];np.savez_compressed(a.out/'priced.npz',H=np.array([v[0] for v in recs],np.int16).reshape(-1,4,24),offset=np.array([v[1] for v in recs],np.int16).reshape(-1,4),parents=np.array([v[2] for v in recs],np.int32).reshape(-1,4),phase=np.array([v[3] for v in recs],np.int16),catalog_index=np.array([v[4] for v in recs]),blocker_indptr=np.r_[0,np.cumsum(list(map(len,bs)))],blocker_indices=np.concatenate(bs) if bs else np.array([],int));report=dict(cycles=n,states=len(kept),histogram=hist,pairs_tried=pairs_tried,seconds=time.monotonic()-start);(a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('FINISHED',json.dumps(report),flush=True)
if __name__=='__main__':main()
