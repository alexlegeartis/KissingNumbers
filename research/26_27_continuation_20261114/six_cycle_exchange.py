"""Exact four-cycle exchange search in the newly enumerated six-support family."""
import argparse,itertools,time,json,sys
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from image_search import above
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'common'))
from golay import golay24

def main():
 ap=argparse.ArgumentParser();ap.add_argument('dim',type=int);ap.add_argument('pool',type=Path);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(exist_ok=False);start=time.monotonic();r=(Path(__file__).resolve().parent/'data');d=a.dim;z=np.load(r/f'PackingStar/25D-31D/new_bounds_20261111/D{d}_{201180 if d==26 else 209862}_coordinates.npz')
 if d==26:O=2*np.r_[z['rational_points'],z['quadratic_A']];OB=2*np.r_[np.zeros_like(z['rational_points']),z['quadratic_B']]
 else:O=2*z['points'];OB=np.zeros_like(O)
 norm=np.sum(O[:,:24].astype(np.int64)**2+3*OB[:,:24].astype(np.int64)**2,axis=1);J=np.setdiff1d(np.arange(24),[0,2,4,5,6,10,11,23]);normal=~np.any(O[:,J],axis=1)&~np.any(OB[:,J],axis=1);fix=(norm!=128)&((norm!=96)|normal);F=O[fix];FB=OB[fix];L=np.load(r/'leech.npy');lk={v.tobytes():i for i,v in enumerate(L)};canon=np.array([i if v[np.flatnonzero(v)[0]]>0 else lk[(-v).tobytes()] for i,v in enumerate(L)]);z=np.load(a.pool);C=z['candidates'];ptr=z['blocker_indptr'];idx=z['blocker_indices'];points={v.tobytes():(v,idx[ptr[i]:ptr[i+1]]) for i,v in enumerate(C)};palette=[]
 for i,j in itertools.combinations(range(d-24),2):
  for s in [-1,1]:
   t=np.zeros(d-24,np.int16);t[[i,j]]=[4,4*s];palette.append(t)
 colorlookup={v.tobytes():i for i,v in enumerate(palette)};bundles={}
 oldids=np.flatnonzero((np.sum(O[:,:24].astype(np.int64)**2,axis=1)==96)&(~np.any(OB,axis=1))&np.all((O[:,:24]==0)|(np.abs(O[:,:24])==4),axis=1))
 for i in oldids:points[O[i].tobytes()]=(O[i],np.array([i],int))
 G=golay24();completions={}
 for row in G[G.sum(1)==8]:
  pos=tuple(map(int,np.flatnonzero(row)))
  for support in itertools.combinations(pos,6):assert support not in completions;completions[support]=pos
 keys=set()
 for c,_ in points.values():
  H=c[:24].copy();T=c[24:].copy()
  if H[np.flatnonzero(H)[0]]<0:H=-H
  if T[np.flatnonzero(T)[0]]<0:T=-T
  keys.add(np.r_[H,T].tobytes())
 edges=[];adj=defaultdict(list)
 for key in keys:
  v=np.frombuffer(key,np.int16);H=v[:24];T=v[24:];P=[np.r_[s*H,t*T].astype(np.int16) for s,t in itertools.product([-1,1],repeat=2)]
  if any(c.tobytes() not in points for c in P):continue
  bs=[points[c.tobytes()][1] for c in P];support=tuple(map(int,np.flatnonzero(H)));octad=completions[support];zeros=np.setdiff1d(octad,support);u=(H//2).copy();u[zeros[0]]=2;u[zeros[1]]=2*int(np.prod(H[list(support)]//4));w=u.copy();w[zeros]*=-1;ia=int(canon[lk[u.tobytes()]]);ib=int(canon[lk[w.tobytes()]]);color=colorlookup[T.tobytes()];e=len(edges);edges.append((ia,ib,color,np.array(P),bs,H));adj[ia].append((ib,e));adj[ib].append((ia,e))
 print('LOADED',len(edges),'edges',len(adj),'parent_lines',flush=True);paths=defaultdict(list);seen=set();cycles=0;best=0;minremoved=10000
 for b,neighbors in adj.items():
  for (va,e),(vc,f) in itertools.combinations(neighbors,2):
   if va==vc or edges[e][2]!=(edges[f][2]^1):continue
   if abs(int(edges[e][5].astype(int)@edges[f][5]))>64:continue
   if va>vc:va,vc=vc,va;e,f=f,e
   key=(va,vc);phase=edges[e][2]
   for bb,ee,ff,pp in paths[key]:
    if bb==b or pp!=(phase^1):continue
    ids=(e,f,ee,ff);signature=tuple(sorted(ids))
    if signature in seen:continue
    seen.add(signature);CA=np.concatenate([edges[j][3] for j in ids]);D=CA.astype(int)@CA.T;np.fill_diagonal(D,-128)
    if D.max()>64:continue
    pp=np.array([va,vc,b,bb]);pc=L[pp];blocked=above(2*pc.astype(np.float32)@F[:,:24].astype(np.float32).T,2*pc.astype(np.float32)@FB[:,:24].astype(np.float32).T).any();badfixed=above(CA.astype(np.float32)@F.astype(np.float32).T,CA.astype(np.float32)@FB.astype(np.float32).T).any()
    if not blocked and not badfixed:
     HH=np.array([edges[j][5] for j in ids]);co=np.array([edges[j][2] for j in ids]);order=np.lexsort(HH.T[::-1]);HH=HH[order];co=co[order];off=co%2^(int(co[0])%2);bk=np.c_[HH,off].astype(np.int16).tobytes();bundles[bk]=(HH,off,pc)
    bs=[x for j in ids for x in edges[j][4]];old=np.unique(np.concatenate(bs));minremoved=min(minremoved,len(old));cycles+=1
    where={int(q):j for j,q in enumerate(old)};M=np.zeros((16,len(old)),bool)
    for j,ns in enumerate(bs):M[j,[where[int(q)] for q in ns]]=True
    graph=csr_matrix(M);matching=maximum_bipartite_matching(graph,perm_type='column');gain=16-int(np.count_nonzero(matching>=0));best=max(best,gain)
    if gain>0:
     reverse={int(q):j for j,q in enumerate(matching) if q>=0};left=set(map(int,np.flatnonzero(matching<0)));right=set();queue=list(left)
     for j in queue:
      for kk in graph.indices[graph.indptr[j]:graph.indptr[j+1]]:
       q=int(kk)
       if matching[j]==q or q in right:continue
       right.add(q)
       if q in reverse and reverse[q] not in left:left.add(reverse[q]);queue.append(reverse[q])
     selected=CA[sorted(left)];removed=old[sorted(right)];keep=np.ones(len(O),bool);keep[removed]=False;RA=np.r_[O[keep],selected];RB=np.r_[OB[keep],np.zeros_like(selected)];assert len(RA)==len(O)+gain
     for j,c in enumerate(selected):
      fail=above(RA.astype(np.int64)@c.astype(np.int64),RB.astype(np.int64)@c.astype(np.int64));fail[len(O)-len(removed)+j]=False;assert not fail.any()
     np.savez_compressed(a.out/'candidate.npz',rational=RA,radical=RB,radicand=3,squared_norm=128,inner_product_limit=64,coordinate_denominator=2);np.savez_compressed(a.out/'exchange.npz',added=selected,removed_indices=removed);print('POSITIVE_EXACT_CHANGED_PAIRS',len(RA),'gain',gain,flush=True);return
   paths[key].append((b,e,f,phase))
  if cycles and cycles%1000==0:print('PROGRESS',cycles,'minremoved',minremoved,'seconds',time.monotonic()-start,flush=True)
 vals=list(bundles.values());np.savez_compressed(a.out/'cycles.npz',H=np.array([v[0] for v in vals],np.int16).reshape(-1,4,24),offset=np.array([v[1] for v in vals],np.int16).reshape(-1,4),parent_vectors=np.array([v[2] for v in vals],np.int16).reshape(-1,4,24));report=dict(edges=len(edges),cycles=cycles,admitted_cycles=len(bundles),min_full_removals=minremoved,best_gain=best,seconds=time.monotonic()-start);(a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('FINISHED',json.dumps(report),flush=True)
if __name__=='__main__':main()
