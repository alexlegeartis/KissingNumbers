"""Enumerate contact squares between Golay octads intersecting in four places."""
import argparse,itertools,time,json,sys
from pathlib import Path
import numpy as np
from image_search import above
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'common'))
from golay import golay24

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--dim',type=int,default=26);ap.add_argument('--pairs',type=int,default=106260);ap.add_argument('--seed',type=int,default=307);a=ap.parse_args();a.out.mkdir(exist_ok=False);r=(Path(__file__).resolve().parent/'data');start=time.monotonic();d=a.dim;z=np.load(r/f'PackingStar/25D-31D/new_bounds_20261111/D{d}_{201180 if d==26 else 209862}_coordinates.npz')
 if d==26:O=2*np.r_[z['rational_points'],z['quadratic_A']];OB=2*np.r_[np.zeros_like(z['rational_points']),z['quadratic_B']]
 else:O=2*z['points'];OB=np.zeros_like(O)
 norm=np.sum(O[:,:24].astype(np.int64)**2+3*OB[:,:24].astype(np.int64)**2,axis=1);J=np.setdiff1d(np.arange(24),[0,2,4,5,6,10,11,23]);normal=~np.any(O[:,J],axis=1)&~np.any(OB[:,J],axis=1);fix=(norm!=128)&((norm!=96)|normal);F=O[fix];FB=OB[fix];FC=F[:,:24].astype(np.float32);BC=FB[:,:24].astype(np.float32);FT=F[:,24:].astype(np.float32);head=(norm==64);hcrows={}
 for v in O[head]:
  h=(v[:24]//2).copy()
  if h[np.flatnonzero(h)[0]]<0:h=-h
  key=h.tobytes();axis=int(np.flatnonzero(v[24:])[0]);hcrows[key]=(h,hcrows.get(key,(None,0))[1]|(1<<axis))
 HC=np.array([v[0] for v in hcrows.values()],np.float32);hflags=np.array([v[1] for v in hcrows.values()]);planepairs=list(itertools.combinations(range(d-24),2));headlimit=np.array([np.where(hflags&((1<<i)|(1<<j)),16,32) for i,j in planepairs]);L=np.load(r/'leech.npy');lk={v.tobytes():i for i,v in enumerate(L)};positive=L[np.arange(len(L)),np.argmax(L!=0,axis=1)]>0;canon=np.array([i if positive[i] else lk[(-v).tobytes()] for i,v in enumerate(L)]);blocked=np.zeros(len(L),bool)
 for at in range(0,len(F),64):blocked|=above(2*L.astype(np.float32)@FC[at:at+64].T,2*L.astype(np.float32)@BC[at:at+64].T).any(1)
 G=golay24();octads=G[G.sum(1)==8];oi,oj=np.where(np.triu(octads.astype(np.int16)@octads.T==4,1));pairs=np.c_[oi,oj];np.random.default_rng(a.seed).shuffle(pairs);Ksign=np.c_[np.ones(8,np.int16),np.array(list(itertools.product([-1,1],repeat=3)),np.int16)];twos=np.array(list(itertools.product([-1,1],repeat=2)),np.int16);outside=np.empty((8,4,4),np.int16);outside[:,:,0]=1;outside[:,:,1:3]=twos[None,:,:];outside[:,:,3]=np.prod(Ksign,axis=1)[:,None]*np.prod(twos,axis=1)[None,:];T=[]
 for i,j in planepairs:
  for s in [-1,1]:
   v=np.zeros(d-24,np.int16);v[[i,j]]=[4,4*s];T.append(v)
 T=np.array(T);offs=np.array([0,1,0,1]);tailcost=np.array([np.abs(T[offs^ph].astype(np.float32)@FT.T) for ph in range(len(T))]);bodyS=np.array([[1,1],[-1,1],[-1,-1],[1,-1]],np.int16);pool={};headpassed=0;tested=0
 for number,(i,j) in enumerate(pairs[:a.pairs]):
  K=np.flatnonzero(octads[i]&octads[j]);A=np.flatnonzero(octads[i]&~octads[j]);B=np.flatnonzero(octads[j]&~octads[i]);assert len(K)==len(A)==len(B)==4;KC=4*Ksign.astype(np.float32)@HC[:,K].T;LA=2*outside.astype(np.float32)@HC[:,A].T;RB=2*outside.astype(np.float32)@HC[:,B].T;worst=np.abs(KC[:,None,None,:])+np.abs(LA[:,:,None,:])+np.abs(RB[:,None,:,:]);good=np.any(np.all(worst[None,:,:,:,:]<=headlimit[:,None,None,None,:],axis=4),axis=0);which=np.argwhere(good);headpassed+=len(which);tested+=128
  if len(which):
   nh=len(which);HH=np.zeros((nh,4,24),np.int16);PP=np.zeros((nh,4,24),np.int16)
   for q,(s,l,t) in enumerate(which):
    HH[q,:,K]=(4*Ksign[s])[:,None];HH[q,:,A]=(2*outside[s,l])[:,None]*bodyS[:,0][None,:];HH[q,:,B]=(2*outside[s,t])[:,None]*bodyS[:,1][None,:];PP[q,:,K]=(2*Ksign[s])[:,None];PP[q,:2,A]=(2*outside[s,l])[:,None]*np.array([1,-1])[None,:];PP[q,2:,B]=(2*outside[s,t])[:,None]*np.array([1,-1])[None,:]
   assert np.all(np.sum(HH.astype(int)**2,axis=2)==96) and np.all(np.sum(PP.astype(int)**2,axis=2)==32)
   pc=np.abs(HH.reshape(-1,24).astype(np.float32)@FC.T).reshape(nh,4,-1);qc=np.abs(HH.reshape(-1,24).astype(np.float32)@BC.T).reshape(nh,4,-1);allow=np.array([~above(pc+tailcost[ph][None,:,:],qc).any((1,2)) for ph in range(len(T))]).T
   for q in np.flatnonzero(allow.any(1)):
    parentids=np.array([int(canon[lk[v.tobytes()]]) for v in PP[q]]);assert len(np.unique(parentids))==4
    if blocked[parentids].any():continue
    H=HH[q];sg=np.where(H[np.arange(4),np.argmax(H!=0,axis=1)]>0,1,-1);H=H*sg[:,None];order=np.lexsort(H.T[::-1]);H=H[order];offset=offs[order];offset=offset^offset[0];key=np.c_[H,offset].astype(np.int16).tobytes();pool[key]=(H,offset,L[parentids])
  if number%500==0:
   report=dict(octad_pairs=number,total_pairs=min(len(pairs),a.pairs),tested_cycles=tested,head_admitted=headpassed,cycles=len(pool),seconds=time.monotonic()-start);print(json.dumps(report),flush=True);(a.out/'progress.json').write_text(json.dumps(report,indent=2)+'\n')
 vals=list(pool.values());np.savez_compressed(a.out/'cycles.npz',H=np.array([v[0] for v in vals],np.int16).reshape(-1,4,24),offset=np.array([v[1] for v in vals],np.int16).reshape(-1,4),parent_vectors=np.array([v[2] for v in vals],np.int16).reshape(-1,4,24));report=dict(octad_pairs=min(len(pairs),a.pairs),tested_cycles=tested,head_admitted=headpassed,cycles=len(pool),seconds=time.monotonic()-start);(a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('FINISHED',json.dumps(report),flush=True)
if __name__=='__main__':main()
