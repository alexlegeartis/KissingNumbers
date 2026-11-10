from adapted_common import check_pattern
from section_replay import validate_claim
"""Independent reconstruction of adaptive sign permutations and selected sections."""
import json,numpy as np
from pathlib import Path
from collections import Counter
import argparse
P=Path(__file__).resolve().parents[1];S=P
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
assert __debug__
z=np.load(S/'data/supports.npz');G=np.array(json.load(open(S/'data/z8_generator.json'))['generator']);coeff=((np.arange(16384)[:,None]>>np.arange(14))&1);C=(coeff@z['sign_generator'])%2
assert len(np.unique(C,axis=0))==16384 and np.min(C[1:].sum(axis=1))==6
patterns=json.load(open(P/'data/adapted_patterns.json'));actual={}
for key,rec in patterns.items():
 check_pattern(key,rec,False)
 a,b=map(int,key.split(','));pos=rec['positions'];flat=sum(pos,[])
 assert len(pos)==a+b and [len(v) for v in pos]==[2]*a+[4]*b and len(set(flat))==len(flat) and all(0<=j<24 for j in flat)
 # Unused sign positions are placed in ascending order into all unselected support cells.
 perm=flat+[j for j in range(24) if j not in flat];assert sorted(perm)==list(range(24))
 ok=np.ones(len(C),dtype=bool)
 for group in pos:ok &=C[:,group].sum(axis=1)*2==len(group)
 count=int(ok.sum());assert count==rec['count'];actual[key]=count
baseline={r['dimension']:r for r in json.load(open(P/'data/section_receipt.json'))}
result=[]
for claim in json.load(open(P/'data/adapted_zero.json')):
 d,sparse,mask=claim['dimension'],claim['terms']['support24'],claim['mask'];T=[i for i in range(24) if mask>>i&1];hist=Counter();total=0
 for lo,hi in zip(z['supports_low'],z['supports_high']):
  support=int(lo)|(int(hi)<<64);ws=[((support>>(4*t))&15).bit_count() for t in T]
  if any(w%2 for w in ws):continue
  key=f'{ws.count(2)},{ws.count(4)}';hist[key]+=1;total+=actual[key]
 assert total==sparse
 image={(0,)*len(T)}
 for grow in G[:,T]:image={tuple((int(x)+a*int(y))%4 for x,y in zip(v,grow)) for v in image for a in range(4)}
 good=sum(all(v) for v in image);assert (8**12*good)%len(image)==0;dense=8**12*good//len(image)
 q=9 if d==88 else 11 if d<92 else 13;assert 6*q<=96-4*len(T)
 terms=dict(dense=dense,transversal24=0,support24=sparse,support6=32*q**4,pure=192+4*len(T));value=sum(terms.values())
 rec=dict(dimension=d,rows=T,mask=mask,terms=terms,total=value,baseline=baseline[d]['total'],gain=value-baseline[d]['total'],projection_size=len(image),nonzero_images=good,pattern_histogram=dict(hist));validate_claim(dict(rec,section='equal' if False else 'zero'),terms);assert rec==claim;result.append(rec);print({'dimension':d,'total':value},flush=True)
args.output.write_text(json.dumps(result,indent=2)+'\n')
