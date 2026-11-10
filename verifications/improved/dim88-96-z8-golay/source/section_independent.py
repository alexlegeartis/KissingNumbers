import json,numpy as np
from pathlib import Path
from collections import Counter
import argparse
from section_replay import finish
P=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
assert __debug__
z=np.load(P/'data/supports.npz'); G=np.array(json.loads((P/'data/z8_generator.json').read_text())['generator'],dtype=np.int64)
# Independent explicit binary matrix multiplication, not imported source proof code.
coeff=((np.arange(1<<14)[:,None]>>np.arange(14))&1).astype(np.int16)
signs=1-2*((coeff@z['sign_generator'])%2)
assert np.unique(signs,axis=0).shape[0]==1<<14
prefix=np.column_stack([np.zeros(len(signs),dtype=np.int16),np.cumsum(signs,axis=1)])
profiles=Counter()
for a,b in zip(z['supports_low'],z['supports_high']):
 s=int(a)+(int(b)<<64); profiles[tuple(((s>>(4*i))&15).bit_count() for i in range(24))]+=1
records=[]
for rec in json.loads((P/'data/sections.json').read_text()):
 if True:
  equal=rec['section']=='equal'
  T=rec['rows']; m=len(T); n=0
  for prof,mult in profiles.items():
   starts=np.cumsum((0,)+prof); sums=np.column_stack([prefix[:,starts[t+1]]-prefix[:,starts[t]] for t in T]); ok=np.all(sums==sums[:,[0]],axis=1) if equal else np.all(sums==0,axis=1)
   n+=mult*int(ok.sum())
  assert n==rec['terms']['support24'],(rec,n)
  # Projection image by integer-coded group closure with unique array rows.
  modulus=8 if equal else 4
  image=np.zeros((1,m),dtype=np.int64)
  for row in G[:,T]:
   image=np.unique(np.concatenate([(image+a*row)%modulus for a in range(modulus)]),axis=0)
  if equal:
   grey=np.array([[(a>>2&1)^((a&1)*x)^((a>>1&1)*y) for x in (0,1) for y in (0,1)] for a in range(8)])
   S=(1-2*grey).sum(axis=1); F=S-2*(1-2*grey[:,0]); count=int(np.all(S[image]==S[image][:,[0]],axis=1).sum())+int(np.all(F[image]==F[image][:,[0]],axis=1).sum())
  else: count=int(np.all(image!=0,axis=1).sum())
  dense=8**12*count//len(image); assert dense==rec['terms']['dense']
  assert 8**12*count%len(image)==0
  assert rec['dimension']==96-m+int(equal)
  assert rec['total']==sum(rec['terms'].values()) and rec['gain']==rec['total']-rec['baseline']
  records.append(finish(rec,dense,n,(1-signs)//2))
args.output.write_text(json.dumps(records,indent=2)+'\n')
print(json.dumps({'status':'PASS','independent_sections':len(records)}))
