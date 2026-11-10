from adapted_common import place_signs
import json,numpy as np,hashlib
from pathlib import Path
from collections import Counter
import argparse
S=Path(__file__).resolve().parents[1];P=S/'data'
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
assert __debug__
z=np.load(S/'data/supports.npz');G=np.array(json.loads((S/'data/z8_generator.json').read_text())['generator'],dtype=np.int64)
labels=np.arange(1<<14);C=((labels[:,None]>>np.arange(14))&1)@z['sign_generator']%2;V=1-2*C
assert len(np.unique(C,axis=0))==16384 and C[1:].sum(axis=1).min()==6
patterns=json.loads((P/'adapted_equal_patterns.json').read_text());counts={}
for key,pat in patterns.items():
 mult=list(map(int,key.split(',')));weights=[w for w,n in enumerate(mult) for _ in range(n)]
 if len(set(w%2 for w in weights))>1:assert pat['count']==0;counts[key]=0;continue
 groups=pat['positions'];flat=sum(groups,[])
 assert list(map(len,groups))==weights and len(set(flat))==len(flat) and all(0<=j<24 for j in flat)
 B=np.zeros((len(groups),24),dtype=np.int64)
 for i,group in enumerate(groups):B[i,group]=1
 sums=V@B.T;count=int(np.all(sums==sums[:,[0]],axis=1).sum());assert count==pat['count'];counts[key]=count
reports=[]
for rec in json.loads((P/'adapted_equal.json').read_text()):
 T=rec['rows'];m=len(T);assert len(set(T))==m and m-1==96-rec['dimension'];hist=Counter();total=0;maps=0;skipped=0
 for lo,hi in zip(z['supports_low'],z['supports_high']):
  mask=int(lo)|(int(hi)<<64);cells=[i for i in range(96) if mask>>i&1];perrow={t:[i for i in cells if i//4==t] for t in T};sizes=[len(perrow[t]) for t in T];key=','.join(str(sizes.count(i)) for i in range(5));hist[key]+=1;total+=counts[key]
  if len(set(w%2 for w in sizes))>1:skipped+=1;continue
  rows=sorted(T,key=lambda t:(len(perrow[t]),t));groups=patterns[key]['positions'];mapping={}
  for t,group in zip(rows,groups):
   assert len(perrow[t])==len(group)
   mapping.update(zip(perrow[t],group))
  mapping.update(zip(sorted(set(cells)-set(mapping)),sorted(set(range(24))-set(mapping.values()))))
  assert mapping==place_signs(cells,T,groups,True)
  assert set(mapping)==set(cells) and sorted(mapping.values())==list(range(24))
  assert [[mapping[c] for c in perrow[t]] for t in rows]==groups
  maps+=1
 assert hist==Counter(rec['pattern_histogram']) and total==rec['terms']['support24']
 image=np.zeros((1,m),dtype=np.int64)
 for row in G[:,T]:image=np.unique(np.concatenate([(image+a*row)%8 for a in range(8)]),axis=0)
 gray=np.array([[1-2*((a>>2&1)^((a&1)*x)^((a>>1&1)*y)) for x in (0,1) for y in (0,1)] for a in range(8)])
 sums=gray.sum(axis=1)[image];flips=(gray.sum(axis=1)-2*gray[:,0])[image];good=int(np.all(sums==sums[:,[0]],axis=1).sum()+np.all(flips==flips[:,[0]],axis=1).sum());assert len(image)==rec['projection_size'] and good==rec['valid_images']
 dense=8**12*good//len(image);assert 8**12*good%len(image)==0
 trans=int(np.all(C[:,T]==C[:,T][:,[0]],axis=1).sum())*4**9
 terms=dict(dense=dense,transversal24=trans,support24=total,support6=32*13**4,pure=192+4*m);assert terms==rec['terms'] and sum(terms.values())==rec['total'] and 78<=96-4*m and rec['total']-rec['baseline']==rec['gain']
 U=np.zeros((m,96),dtype=np.int64)
 for i,t in enumerate(T):U[i,4*t:4*t+4]=1
 A=U[1:]-U[0];assert np.array_equal(A[:,[4*t for t in T[1:]]],np.eye(m-1,dtype=np.int64))
 out=dict(dimension=rec['dimension'],total=rec['total'],baseline=rec['baseline'],gain=rec['gain'],terms=terms,physical_permutations=maps,mixed_parity_supports=skipped,projection_size=len(image),valid_images=good,section_rank=m-1);reports.append(out);print(out,flush=True)
args.output.write_text(json.dumps({'status':'PASS','dimensions':reports},indent=2)+'\n')
