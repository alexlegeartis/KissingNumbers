from guards import validate_source,validate_pattern
import json,numpy as np,hashlib,itertools
from pathlib import Path
from collections import Counter
import argparse
ROOT=Path(__file__).resolve().parents[1];P=ROOT/'data';S=ROOT.parent
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
assert __debug__
data=json.loads((P/'nr_partition.json').read_text());G=data['octacode_generator'];assert len(G)==4 and all(len(row)==8 and all(isinstance(x,int) and 0<=x<4 for x in row) for row in G)
gray=[0,1,3,2]
NR={sum(gray[sum(c[i]*G[i][j] for i in range(4))%4]<<(2*j) for j in range(8)) for c in itertools.product(range(4),repeat=4)}
assert len(NR)==256 and NR==set(data['NR_words']) and len(data['NR_words'])==256
nrl=list(NR);nrhist=Counter((x^y).bit_count() for i,x in enumerate(nrl) for y in nrl[:i]);assert min(nrhist)==6
reps=data['translate_representatives'];assert len(reps)==8 and all(0<=t<2**16 for t in reps)
parts=[{x^t for x in NR} for t in reps];assert all(parts[i]==set(data['partition'][i]) for i in range(8))
H=set.union(*parts);assert len(H)==2048 and sum(map(len,parts))==2048
# Independently reconstruct claimed linear span, rather than assuming union is linear.
span={0}
for row in data['span_basis']:span |= {x^row for x in span}
assert len(span)==2048 and span==H and NR<=span and min(x.bit_count() for x in span if x)==4
W=[]
for label in range(16):
 color=(label&7)^(7 if label&8 else 0);head=sum((1 if label>>j&1 else 2)<<(2*j) for j in range(4));W +=[head|(w<<8) for w in parts[color]]
source=json.loads((P/'nr_pair4.json').read_text());assert len(W)==len(set(W))==4096 and set(W)==set(source) and len(source)==4096
validate_source(data,source,NR,parts,W)
hist=Counter((x^y).bit_count() for i,x in enumerate(W) for y in W[:i]);assert sum(hist.values())==8386560 and min(hist)==6
print('NR construction, partition, span and scalar all-pairs PASS',dict(hist),flush=True)
z=np.load(S/'data/supports.npz');linear=[0]
for row in z['sign_generator']:
 mask=sum(int(b)<<j for j,b in enumerate(row));linear +=[w^mask for w in linear]
assert len(set(linear))==16384 and min(w.bit_count() for w in linear if w)==6
R=np.array(json.loads((S/'data/z8_generator.json').read_text())['generator'],dtype=np.int64)
reports=[]
for equal,file in [(False,'nonlinear_patterns.json'),(True,'equal_nonlinear_patterns.json')]:
 patterns=json.loads((P/file).read_text());counts={}
 for key,pat in patterns.items():
  validate_pattern(key,pat,equal)
  sizes=list(map(int,key.split(',')));weights=[i for i,n in enumerate(sizes) for _ in range(n)] if equal else [2]*sizes[0]+[4]*sizes[1]
  assert pat['family'] in ('NR','linear')
  if equal and len(set(w%2 for w in weights))>1:assert pat['count']==0;counts[key]=0;continue
  groups=pat['positions'];flat=sum(groups,[]);assert list(map(len,groups))==weights and len(set(flat))==len(flat) and all(0<=j<24 for j in flat)
  translation=pat.get('translation',0);assert 0<=translation<2**24 and (pat['family']=='NR' or translation==0)
  masks=[sum(1<<j for j in group) for group in groups];words=linear if pat['family']=='linear' else W
  def good(w):
   sums=[len(group)-2*((w^translation)&mask).bit_count() for group,mask in zip(groups,masks)]
   return len(set(sums))==1 if equal else not any(sums)
  n=sum(good(w) for w in words);assert n==pat['count'];counts[key]=n
 for rec in json.loads((P/'verified.json').read_text())['records']:
  if (rec['section']=='equal')!=equal:continue
  T=rec['rows'];m=len(T);assert len(set(T))==m;histp=Counter();families=Counter();total=maps=0
  for lo,hi in zip(z['supports_low'],z['supports_high']):
   support=int(lo)|(int(hi)<<64);cells=[j for j in range(96) if support>>j&1];byrow={t:[j for j in cells if j//4==t] for t in T};ws=[len(byrow[t]) for t in T]
   if not equal and any(w%2 for w in ws):continue
   key=','.join(str(ws.count(i)) for i in range(5)) if equal else f'{ws.count(2)},{ws.count(4)}';pat=patterns[key];histp[key]+=1;total+=counts[key];families[pat['family']]+=counts[key]
   if equal and len(set(w%2 for w in ws))>1:continue
   rows=sorted((t for t in T if equal or byrow[t]),key=lambda t:(len(byrow[t]),t));mapping={}
   for t,g in zip(rows,pat['positions']):
    assert len(byrow[t])==len(g)
    mapping.update(zip(byrow[t],g))
   mapping.update(zip(sorted(set(cells)-set(mapping)),sorted(set(range(24))-set(mapping.values()))))
   assert set(mapping)==set(cells) and sorted(mapping.values())==list(range(24)) and [[mapping[c] for c in byrow[t]] for t in rows]==pat['positions'];maps+=1
  assert histp==Counter(rec['pattern_histogram']) and families==Counter(rec['family_points']) and total==rec['terms']['support24'] and maps==rec['physical_permutations']
  mod=8 if equal else 4;image=np.zeros((1,m),dtype=np.int64)
  for row in R[:,T]:image=np.unique(np.concatenate([(image+a*row)%mod for a in range(mod)]),axis=0)
  if equal:
   gray8=np.array([[1-2*((a>>2&1)^((a&1)*x)^((a>>1&1)*y)) for x in (0,1) for y in (0,1)] for a in range(8)])
   sums=gray8.sum(axis=1)[image];flips=(gray8.sum(axis=1)-2*gray8[:,0])[image];good=int(np.all(sums==sums[:,[0]],axis=1).sum()+np.all(flips==flips[:,[0]],axis=1).sum())
  else:good=int(np.all(image!=0,axis=1).sum())
  assert good==rec['valid_images'] and len(image)==rec['projection_size'] and 8**12*good%len(image)==0
  dense=8**12*good//len(image);trans=sum(len({w>>t&1 for t in T})==1 for w in linear)*4**9 if equal else 0
  d=rec['dimension'];q=9 if d==88 else 11 if d<=91 else 13;assert 6*q<=96-4*m
  U=np.zeros((m,96),dtype=np.int64)
  for i,t in enumerate(T):U[i,4*t:4*t+4]=1
  A=U[1:]-U[0] if equal else U;assert len(A)==96-d and np.array_equal(A[:,[4*t for t in (T[1:] if equal else T)]],np.eye(len(A),dtype=np.int64))
  terms=dict(dense=dense,transversal24=trans,support24=total,support6=32*q**4,pure=192+4*m);assert terms==rec['terms'] and sum(terms.values())==rec['total'] and rec['total']-rec['baseline']==rec['gain'];reports.append(rec);print(d,rec['total'],rec['gain'],'PASS',flush=True)
args.output.write_text(json.dumps(dict(status='PASS',NR24_pairs=8386560,records=reports),indent=2)+'\n')
