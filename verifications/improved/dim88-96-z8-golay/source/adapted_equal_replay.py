from adapted_common import check_pattern
from section_replay import validate_claim
import json,numpy as np
from pathlib import Path
from collections import Counter
import argparse
P=Path(__file__).resolve().parents[1];S=P
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
assert __debug__
z=np.load(S/'data/supports.npz');G=np.array(json.load(open(S/'data/z8_generator.json'))['generator'])
# Build words by XOR doubling, independent of the pattern search matrix expansion.
W=[0]
for row in z['sign_generator']:
 v=sum(int(b)<<i for i,b in enumerate(row));W +=[w^v for w in W]
actual={}
for key,rec in json.load(open(P/'data/adapted_equal_patterns.json')).items():
 check_pattern(key,rec,True)
 sizes=[int(a) for a in key.split(',')];weights=sum(([i]*n for i,n in enumerate(sizes)),[])
 if (sizes[1]+sizes[3]) and (sizes[0]+sizes[2]+sizes[4]):assert rec['count']==0;actual[key]=0;continue
 groups=rec['positions'];flat=sum(groups,[]);assert [len(g) for g in groups]==weights and len(flat)==len(set(flat)) and all(0<=j<24 for j in flat)
 masks=[sum(1<<j for j in g) for g in groups]
 count=sum(len({len(g)-2*(w&mask).bit_count() for g,mask in zip(groups,masks)})==1 for w in W)
 assert count==rec['count'];actual[key]=count
base={r['dimension']:r for r in json.load(open(P/'data/section_receipt.json'))};result=[]
for claim in json.load(open(P/'data/adapted_equal.json')):
 d,sparse,T=claim['dimension'],claim['terms']['support24'],claim['rows'];hist=Counter();total=0
 for lo,hi in zip(z['supports_low'],z['supports_high']):
  support=int(lo)|(int(hi)<<64);ws=[((support>>(4*t))&15).bit_count() for t in T];key=','.join(str(ws.count(i)) for i in range(5));hist[key]+=1;total+=actual[key]
 assert total==sparse
 image={(0,)*len(T)}
 for grow in G[:,T]:image={tuple((int(x)+a*int(y))%8 for x,y in zip(v,grow)) for v in image for a in range(8)}
 gray=[[1-2*((a>>2&1)^((a&1)*x)^((a>>1&1)*y)) for x in (0,1) for y in (0,1)] for a in range(8)];sums=[sum(r) for r in gray];flips=[sum(r)-2*r[0] for r in gray]
 good=sum(len({sums[a] for a in v})==1 for v in image)+sum(len({flips[a] for a in v})==1 for v in image);assert 8**12*good%len(image)==0;dense=8**12*good//len(image)
 trans=sum(len({w>>t&1 for t in T})==1 for w in W)*4**9
 q=13;assert 6*q<=96-4*len(T)
 terms=dict(dense=dense,transversal24=trans,support24=sparse,support6=32*q**4,pure=192+4*len(T));value=sum(terms.values());rec=dict(dimension=d,rows=T,terms=terms,total=value,baseline=base[d]['total'],gain=value-base[d]['total'],projection_size=len(image),valid_images=good,pattern_histogram=dict(hist));validate_claim(dict(rec,section='equal' if True else 'zero'),terms);assert rec==claim;result.append(rec);print({'dimension':d,'total':value},flush=True)
args.output.write_text(json.dumps(result,indent=2)+'\n')
