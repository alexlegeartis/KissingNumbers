"""Direct contact-cycle exchanges, restoring parents released by partial deletions."""
import argparse,itertools,time,json
from pathlib import Path
import numpy as np
from image_search import above

def main():
 ap=argparse.ArgumentParser();ap.add_argument('pack',type=Path);ap.add_argument('catalog',type=Path);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--skip',type=int,default=0);a=ap.parse_args();a.out.mkdir(exist_ok=False);start=time.monotonic();r=(Path(__file__).resolve().parent/'data');z=np.load(a.pack/'data.npz');d=int(z['dimension']);T=z['tail_palette'];nt=len(T);L=z['leech'];F=z['fixed_A'];FB=z['fixed_B'];flat=z['state'][z['initial']];ci=flat//nt;phase=flat%nt;OH=z['H'].reshape(-1,4,24)[ci];OC=z['offset'].reshape(-1,4)[ci]^phase[:,None];OP=z['parents'][ci];nold=len(ci);inc=(np.abs(np.einsum('cbi,cpi->cbp',OH.astype(np.int64),L[OP].astype(np.int64)))==48).astype(np.uint8);assert np.all(inc.sum(1)==2) and np.all(inc.sum(2)==2);owner=np.full(len(L),-1,np.int32);owner[OP.ravel()]=np.repeat(np.arange(nold),4);cost=np.abs(T.astype(np.int64)@T.T);oldHF=OH.reshape(-1,24).astype(np.float32);lk={v.tobytes():i for i,v in enumerate(L)};c=np.load(a.catalog);H=c['H'];off=c['offset'];PV=c['parent_vectors'];P=np.array([[lk[(v if v[np.flatnonzero(v)[0]]>0 else -v).tobytes()] for v in ps] for ps in PV],np.int32);n=len(H);best=-100000;mincost=100000;print('LOADED',d,n,'cycles',flush=True)
 def positive(h,colors,parents,badbody,restore,gain):
  baseline=201180 if d==26 else 209862;src=np.load(r/f'PackingStar/25D-31D/new_bounds_20261111/D{d}_{baseline}_coordinates.npz')
  if d==26:O=2*np.r_[src['rational_points'],src['quadratic_A']];OB=2*np.r_[np.zeros_like(src['rational_points']),src['quadratic_B']]
  else:O=2*src['points'];OB=np.zeros_like(O)
  lookup={v.tobytes():j for j,(v,b) in enumerate(zip(O,OB)) if not b.any()};removed=[]
  for cyc,body in np.argwhere(badbody):
   hh=OH[cyc,body];tt=T[OC[cyc,body]]
   for s,t in itertools.product([-1,1],repeat=2):removed.append(lookup[np.r_[s*hh,t*tt].astype(np.int16).tobytes()])
  for p in parents:
   for s in [-1,1]:
    key=np.r_[2*s*L[p],np.zeros(d-24,np.int16)].astype(np.int16).tobytes()
    if key in lookup:removed.append(lookup[key])
  C=np.array([np.r_[s*hh,t*T[int(cc)]] for hh,cc in zip(h,colors) for s,t in itertools.product([-1,1],repeat=2)],np.int16);rest=np.array([np.r_[2*s*L[p],np.zeros(d-24,np.int16)] for p in restore for s in [-1,1]],np.int16).reshape(-1,d);C=np.r_[C,rest];removed=np.unique(removed);keep=np.ones(len(O),bool);keep[removed]=False;RA=np.r_[O[keep],C];RB=np.r_[OB[keep],np.zeros_like(C)];assert len(RA)==len(O)+gain and gain>0;assert len({np.r_[v,b].tobytes() for v,b in zip(RA,RB)})==len(RA);assert np.all(np.sum(RA.astype(np.int64)**2+3*RB.astype(np.int64)**2,axis=1)==128) and not np.any(np.sum(RA.astype(np.int64)*RB.astype(np.int64),axis=1));A64=RA.astype(np.int64);B64=RB.astype(np.int64)
  for j,v in enumerate(C):
   fail=above(A64@v.astype(np.int64),B64@v.astype(np.int64));fail[len(O)-len(removed)+j]=False;assert not fail.any(),('changed geometry',j,np.flatnonzero(fail)[:8])
  np.savez_compressed(a.out/'candidate.npz',rational=RA,radical=RB,radicand=3,squared_norm=128,inner_product_limit=64,coordinate_denominator=2);np.savez_compressed(a.out/'exchange.npz',added_A=C,added_B=np.zeros_like(C),removed_indices=removed,restored_parent_lines=restore,new_H=h,new_colors=colors,new_parents=parents);report=dict(dimension=d,count=len(RA),gain=gain,added=len(C),removed=len(removed),restored_parent_lines=len(restore),seconds=time.monotonic()-start);(a.out/'positive.json').write_text(json.dumps(report,indent=2)+'\n');print('POSITIVE_EXACT_CHANGED_PAIRS',json.dumps(report),flush=True)
 incmask=np.sum(inc*np.left_shift(np.uint8(1),np.arange(4,dtype=np.uint8))[None,:,None],axis=1,dtype=np.uint8)
 for at in range(a.skip,n,128):
  hh=H[at:at+128];oo=off[at:at+128];pp=P[at:at+128];nb=len(hh);D=np.abs(hh.reshape(-1,24).astype(np.float32)@oldHF.T).reshape(nb,4,nold,4);pc=np.abs(hh.reshape(-1,24).astype(np.float32)@F[:,:24].astype(np.float32).T).reshape(nb,4,-1);qc=np.abs(hh.reshape(-1,24).astype(np.float32)@FB[:,:24].astype(np.float32).T).reshape(nb,4,-1);retcount=np.sum(owner[pp]<0,axis=1);covered=np.any(OP[None,:,:,None]==pp[:,None,None,:],axis=3)
  parity=oo[:,:,None,None]^(OC%2).reshape(1,1,nold,4);weak=D>32;weak0=(weak&(parity==0)).any(1);weak1=(weak&(parity==1)).any(1);strong=(D>48).any(1);hard=(D>64).any(1)
  for ph in range(nt):
   colors=oo^ph;tc=np.abs(T[colors].astype(np.float32)@F[:,24:].astype(np.float32).T);allow=~above(pc+tc,qc).any((1,2));badbody=np.where((OC//2)[None,:,:]==ph//2,hard|(weak1 if ph&1 else weak0),strong);midcount=badbody.sum((1,2));bm=np.packbits(badbody,axis=2,bitorder='little')[:,:,0];released=((bm[:,:,None]&incmask[None,:,:])==incmask[None,:,:])&~covered;restorecount=released.sum((1,2));removed=4*midcount+2*retcount;gain=16-removed+2*restorecount
   if allow.any():
    vv=int(gain[allow].max());mincost=min(mincost,int(removed[allow].min()))
    if vv>best:best=vv;print('BEST',best,'at',at,'phase',ph,'seconds',time.monotonic()-start,flush=True)
   for j in np.flatnonzero(allow&(gain>0)):
    restore=OP[released[j]];positive(hh[j],colors[j],pp[j],badbody[j],restore,int(gain[j]));return
  if at%4096==0:print('PROGRESS',at,n,'best',best,'minremoved',mincost,'seconds',round(time.monotonic()-start,1),flush=True)
 report=dict(dimension=d,cycles=n,start_index=a.skip,best_gain=best,min_full_removals=mincost,seconds=time.monotonic()-start);(a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('FINISHED',json.dumps(report),flush=True)
if __name__=='__main__':main()
