"""Search connected incumbent-cycle cavities by exact compatible-clique recursion."""
import argparse,itertools,time,json,sys
from pathlib import Path
from collections import defaultdict,deque
import numpy as np
from cycle_reconstruct import reconstruct

def main():
 ap=argparse.ArgumentParser();ap.add_argument('pack',type=Path);ap.add_argument('priced',type=Path,nargs='+');ap.add_argument('--out',type=Path,required=True);ap.add_argument('--max-free',type=int,default=4);ap.add_argument('--seconds',type=float,default=1200);ap.add_argument('--per-scope',type=float,default=.3);a=ap.parse_args();a.out.mkdir(exist_ok=False);start=time.monotonic();z=np.load(a.pack/'data.npz');T=z['tail_palette'];nt=len(T);flat=z['state'][z['initial']];ci=flat//nt;ph=flat%nt;OH=z['H'].reshape(-1,4,24)[ci];OC=z['offset'].reshape(-1,4)[ci]^ph[:,None];OP=z['parents'][ci];nold=len(ci);vals=[];known={}
 def put(h,c,p,b):
  h=h.copy();sg=np.where(h[np.arange(4),np.argmax(h!=0,axis=1)]>0,1,-1);h*=sg[:,None];order=np.lexsort(h.T[::-1]);h=h[order];c=c[order];key=np.c_[h,c].astype(np.int16).tobytes();bb=tuple(sorted(map(int,b)))
  if key in known:assert vals[known[key]][3]==bb;return
  known[key]=len(vals);vals.append((h,c,np.sort(p),bb))
 for i in range(nold):put(OH[i],OC[i],OP[i],[i])
 for path in a.priced:
  z=np.load(path);ptr=z['blocker_indptr'];idx=z['blocker_indices']
  for j,h in enumerate(z['H']):
   bs=idx[ptr[j]:ptr[j+1]]
   if len(bs)<=a.max_free:put(h,z['offset'][j]^int(z['phase'][j]),z['parents'][j],bs)
 H=np.array([v[0] for v in vals],np.int16);C=np.array([v[1] for v in vals],np.int16);P=np.array([v[2] for v in vals],np.int32);B=[v[3] for v in vals];n=len(H);HF=H.reshape(-1,24).astype(np.float32);cost=np.abs(T.astype(np.int64)@T.T);adj=[];print('LOADED',n,'states',nold,'old',flush=True)
 for at in range(0,n,32):
  hh=H[at:at+32];co=C[at:at+32];pp=P[at:at+32];nb=len(hh);D=np.abs(hh.reshape(-1,24).astype(np.float32)@HF.T).reshape(nb,4,n,4);bad=(D+cost[co[:,:,None,None],C.reshape(1,1,n,4)]>64).any((1,3));bad|=np.any(pp[:,None,:,None]==P[None,:,None,:],axis=(2,3));ok=~bad;ok[np.arange(nb),at+np.arange(nb)]=False;packed=np.packbits(ok,axis=1,bitorder='little');adj.extend(int.from_bytes(row.tobytes(),'little') for row in packed)
  if at%2048==0:print('GRAPH',at,n,round(time.monotonic()-start,1),flush=True)
 for i in range(nold):assert (adj[i]&((1<<nold)-1)).bit_count()==nold-1
 bucket=defaultdict(int)
 for j,b in enumerate(B):bucket[b]|=1<<j
 edges=[b for b in bucket if len(b)>=2];incident=defaultdict(list)
 for b in edges:
  for v in b:incident[v].append(b)
 scopes={tuple([v]) for v in range(nold)}|set(edges);queue=deque(edges)
 while queue:
  ss=queue.popleft();near=set(b for v in ss for b in incident[v])
  for b in near:
   u=tuple(sorted(set(ss)|set(b)))
   if len(u)<=a.max_free and u not in scopes:scopes.add(u);queue.append(u)
 def eligible(ss):
  bits=0
  for k in range(1,len(ss)+1):
   for b in itertools.combinations(ss,k):bits|=bucket.get(b,0)
  return bits
 ordered=[]
 for ss in scopes:
  bits=eligible(ss)
  if bits.bit_count()>len(ss):ordered.append((ss,bits))
 ordered.sort(key=lambda v:(len(v[0]),-v[1].bit_count()));print('SCOPES',len(ordered),flush=True);nodes=0;timeouts=0
 def colors_enough(bits,need):
  colors=0;left=bits
  while left:
   colors+=1
   if colors>=need:return True
   available=left
   while available:
    bit=available&-available;v=bit.bit_length()-1;left^=bit;available^=bit;available&=~adj[v]
  return False
 def dfs(bits,need,chosen,deadline):
  nonlocal nodes
  nodes+=1
  if nodes%128==0 and time.monotonic()>deadline:raise TimeoutError
  if not need:return chosen
  if bits.bit_count()<need:return None
  if need==1:return chosen+[(bits&-bits).bit_length()-1]
  if not colors_enough(bits,need):return None
  vertices=[];tmp=bits
  while tmp:
   bit=tmp&-tmp;tmp^=bit;v=bit.bit_length()-1;vertices.append(v)
  vertices.sort(key=lambda v:(v>=nold,(bits&adj[v]).bit_count()),reverse=True)
  for v in vertices:
   bits&=~(1<<v)
   if bits.bit_count()<need-1:return None
   found=dfs(bits&adj[v],need-1,chosen+[v],deadline)
   if found is not None:return found
  return None
 for tr,(ss,bits) in enumerate(ordered):
  if time.monotonic()-start>a.seconds:break
  try:found=dfs(bits,len(ss)+1,[],min(start+a.seconds,time.monotonic()+a.per_scope))
  except TimeoutError:timeouts+=1;found=None
  if found is not None:
   report=reconstruct(a.pack,a.out,H[found],C[found],P[found],ss);np.savez_compressed(a.out/'trade.npz',H=H[found],colors=C[found],parents=P[found],deleted=ss,states=found);print('POSITIVE_EXACT_CHANGED_PAIRS',json.dumps(report),flush=True);return
  if tr%100==0:print('SEARCH',tr,len(ordered),'nodes',nodes,'timeouts',timeouts,'seconds',round(time.monotonic()-start,1),flush=True)
 report=dict(states=n,old=nold,scopes=len(ordered),searched=tr+1 if ordered else 0,nodes=nodes,timeouts=timeouts,seconds=time.monotonic()-start);(a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('FINISHED',json.dumps(report),flush=True)
if __name__=='__main__':main()
