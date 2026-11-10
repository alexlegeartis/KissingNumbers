from adapted_common import place_signs
import json,numpy as np,hashlib
from pathlib import Path
from collections import Counter
import argparse
S=Path(__file__).resolve().parents[1];P=S/'data'
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
assert __debug__
z=np.load(S/'data/supports.npz');G=np.array(json.loads((S/'data/z8_generator.json').read_text())['generator'],dtype=np.int64)
# Independent scalar-mask binary enumeration.
words=[0]
for row in z['sign_generator']:
 mask=sum(int(v)<<j for j,v in enumerate(row));words +=[x^mask for x in words]
assert len(set(words))==16384 and min(x.bit_count() for x in words if x)==6
patterns=json.loads((P/'adapted_patterns.json').read_text());counts={}
for key,pat in patterns.items():
 a,b=map(int,key.split(','));groups=pat['positions'];flat=sum(groups,[])
 assert list(map(len,groups))==[2]*a+[4]*b and len(set(flat))==len(flat) and all(0<=i<24 for i in flat)
 masks=[sum(1<<i for i in group) for group in groups]
 count=sum(all((w&mask).bit_count()*2==mask.bit_count() for mask in masks) for w in words)
 assert count==pat['count'];counts[key]=count
reports=[]
for rec in json.loads((P/'adapted_zero.json').read_text()):
 T=rec['rows'];m=len(T);assert m==96-rec['dimension'] and len(set(T))==m and sum(1<<i for i in T)==rec['mask']
 total=0;hist=Counter();mapped=0;skipped=0
 for lo,hi in zip(z['supports_low'],z['supports_high']):
  support=int(lo)|(int(hi)<<64);cells=[i for i in range(96) if support>>i&1]
  perrow={t:[i for i in cells if i//4==t] for t in T}
  if any(len(v)%2 for v in perrow.values()):skipped+=1;continue
  # Canonical explicit recipe: two-cell selected rows first, then four-cell;
  # ties by physical row index, cells within each row ascending.
  rows=sorted((t for t in T if perrow[t]),key=lambda t:(len(perrow[t]),t))
  a=sum(len(perrow[t])==2 for t in rows);b=sum(len(perrow[t])==4 for t in rows);key=f'{a},{b}';pat=patterns[key]
  mapping={}
  for t,group in zip(rows,pat['positions']):
   for physical,signpos in zip(perrow[t],group):mapping[physical]=signpos
  remaining_positions=sorted(set(range(24))-set(mapping.values()))
  remaining_cells=sorted(set(cells)-set(mapping))
  mapping.update(zip(remaining_cells,remaining_positions))
  assert mapping==place_signs(cells,T,pat["positions"],False)
  assert set(mapping)==set(cells) and sorted(mapping.values())==list(range(24))
  actual_groups=[sum(1<<mapping[c] for c in perrow[t]) for t in rows]
  assert actual_groups==[sum(1<<j for j in group) for group in pat['positions']]
  # Thus every sign pattern accepted above has physical row sums zero,
  # and bijective coordinate placement preserves sign-code distance >=6.
  hist[key]+=1;total+=counts[key];mapped+=1
 assert dict(hist)==rec['pattern_histogram'] and total==rec['terms']['support24']
 image=np.zeros((1,m),dtype=np.int64)
 for row in G[:,T]:image=np.unique(np.concatenate([(image+a*row)%4 for a in range(4)]),axis=0)
 good=int(np.all(image!=0,axis=1).sum());assert len(image)==rec['projection_size'] and good==rec['nonzero_images'];dense=8**12*good//len(image);assert dense==rec['terms']['dense'] and 8**12*good%len(image)==0
 q=9 if rec['dimension']==88 else 13 if rec['dimension']==92 else 11
 terms=dict(dense=dense,transversal24=0,support24=total,support6=32*q**4,pure=192+4*m)
 assert 6*q<=96-4*m and terms==rec['terms'] and sum(terms.values())==rec['total'] and rec['total']-rec['baseline']==rec['gain']
 reports.append(dict(dimension=rec['dimension'],total=rec['total'],baseline=rec['baseline'],gain=rec['gain'],support24=total,actual_physical_permutations=mapped,odd_row_supports_excluded=skipped,projection_size=len(image),good_images=good,section_rank=m))
 print(reports[-1],flush=True)
args.output.write_text(json.dumps({'status':'PASS','dimensions':reports},indent=2)+'\n')
