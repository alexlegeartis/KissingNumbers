from guards import validate_source,validate_pattern
"""Independent NR construction, all pairs, all patterns and physical maps."""
import json,numpy as np
from pathlib import Path
from collections import Counter
import argparse
ROOT=Path(__file__).resolve().parents[1];P=ROOT/'data';S=ROOT.parent;A=S/'data'
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
assert __debug__
z=np.load(S/'data/supports.npz');R=np.array(json.load(open(S/'data/z8_generator.json'))['generator']);data=json.load(open(P/'nr_partition.json'));G=np.array(data['octacode_generator']);coef=np.array([[(i//4**j)%4 for j in range(4)] for i in range(256)]);C=coef@G%4;gray=np.array([[0,0],[1,0],[1,1],[0,1]],dtype=np.int16);NR=gray[C].reshape(256,16);nrwords=[sum(int(v)<<j for j,v in enumerate(row)) for row in NR];assert set(nrwords)==set(data['NR_words']) and len(set(nrwords))==256
nrweights=NR.sum(axis=1);dist=nrweights[:,None]+nrweights[None,:]-2*(NR@NR.T);np.fill_diagonal(dist,100);assert dist.min()==6
parts=[{x^t for x in nrwords} for t in data['translate_representatives']];assert len(parts)==8 and sum(len(x) for x in parts)==len(set.union(*parts))==2048
H=np.array([[(x>>j)&1 for j in range(16)] for x in set.union(*parts)],dtype=np.int16);hw=H.sum(axis=1);dist=hw[:,None]+hw[None,:]-2*(H@H.T);np.fill_diagonal(dist,100);assert dist.min()==4
rebuilt=[]
for label in range(16):
 color=(label&7)^(7 if label&8 else 0);head=sum((1 if label>>j&1 else 2)<<(2*j) for j in range(4));rebuilt +=[head|(w<<8) for w in parts[color]]
source=json.load(open(P/'nr_pair4.json'));assert len(source)==len(set(source))==4096 and set(source)==set(rebuilt)
validate_source(data,source,nrwords,parts,rebuilt)
NR24=np.array([[(x>>j)&1 for j in range(24)] for x in source],dtype=np.int16);weights=NR24.sum(axis=1);minimum=24;pairs=0
for start in range(0,4096,128):
 end=min(start+128,4096);dist=weights[start:end,None]+weights[None,:end]-2*(NR24[start:end]@NR24[:end].T)
 for i in range(end-start):
  j=start+i;pairs+=j
  if j:minimum=min(minimum,int(dist[i,:j].min()))
assert pairs==4096*4095//2 and minimum==6
labels=(np.arange(16384)[:,None]>>np.arange(14))&1;linear=(labels@z['sign_generator'])%2;assert len(np.unique(linear,axis=0))==16384 and linear[1:].sum(axis=1).min()==6
print('Independent NR256/partition/H2048/NR24: PASS',pairs,'NR24 pairs',flush=True)
base={r['dimension']:r for n in ['adapted_zero.json','adapted_equal.json'] for r in json.load(open(A/n))};reports=[]
for equal,patternfile,rowfile in [(False,'nonlinear_patterns.json','final_zero_rows.txt'),(True,'equal_nonlinear_patterns.json','equal_results.txt')]:
 patterns=json.load(open(P/patternfile));counts={}
 for key,pat in patterns.items():
  validate_pattern(key,pat,equal)
  sizes=list(map(int,key.split(',')));wanted=sum(([i]*n for i,n in enumerate(sizes)),[]) if equal else [2]*sizes[0]+[4]*sizes[1]
  if equal and len({w%2 for w in wanted})>1:assert pat['count']==0;counts[key]=0;continue
  groups=pat['positions'];flat=sum(groups,[]);assert list(map(len,groups))==wanted and len(flat)==len(set(flat)) and all(0<=v<24 for v in flat)
  assert pat['family'] in ['linear','NR'];B=linear if pat['family']=='linear' else NR24^np.array([(pat['translation']>>j)&1 for j in range(24)],dtype=np.int16)
  assert 0<=pat.get('translation',0)<2**24
  sums=np.stack([(1-2*B[:,g]).sum(axis=1) for g in groups],axis=1) if groups else np.zeros((len(B),0),dtype=np.int16)
  ok=np.all(sums==sums[:,[0]],axis=1) if equal else np.all(sums==0,axis=1);n=int(ok.sum());assert n==pat['count'];counts[key]=n
 for line in (P/rowfile).read_text().splitlines():
  d,expected,mask=map(int,line.split());T=[j for j in range(24) if mask>>j&1];hist=Counter();total=maps=0;family_points=Counter()
  for lo,hi in zip(z['supports_low'],z['supports_high']):
   support=int(lo)|(int(hi)<<64);cells=[j for j in range(96) if support>>j&1];byrow={t:[j for j in cells if j//4==t] for t in T};ws=[len(byrow[t]) for t in T]
   if not equal and any(w%2 for w in ws):continue
   key=','.join(str(ws.count(i)) for i in range(5)) if equal else f'{ws.count(2)},{ws.count(4)}';pat=patterns[key];hist[key]+=1;total+=counts[key];family_points[pat['family']]+=counts[key]
   if equal and len({w%2 for w in ws})>1:continue
   rows=sorted((t for t in T if equal or byrow[t]),key=lambda t:(len(byrow[t]),t));mapping={}
   for t,g in zip(rows,pat['positions']):mapping.update(zip(byrow[t],g))
   mapping.update(zip(sorted(set(cells)-set(mapping)),sorted(set(range(24))-set(mapping.values()))));assert set(mapping)==set(cells) and sorted(mapping.values())==list(range(24));assert [[mapping[c] for c in byrow[t]] for t in rows]==pat['positions'];maps+=1
  assert total==expected
  mod=8 if equal else 4;image={(0,)*len(T)}
  for grow in R[:,T]:image={tuple((int(x)+a*int(y))%mod for x,y in zip(v,grow)) for v in image for a in range(mod)}
  if equal:
   rowsums=[sum(1-2*((a>>2&1)^((a&1)*x)^((a>>1&1)*y)) for x in (0,1) for y in (0,1)) for a in range(8)];flip=[rowsums[a]-2*(1-2*(a>>2&1)) for a in range(8)];good=sum(len({rowsums[x] for x in v})==1 for v in image)+sum(len({flip[x] for x in v})==1 for v in image)
  else:good=sum(all(v) for v in image)
  assert (8**12*good)%len(image)==0;dense=8**12*good//len(image);trans=int(np.all(linear[:,T]==linear[:,[T[0]]],axis=1).sum())*4**9 if equal else 0;q=9 if d==88 else 11 if d<=91 else 13;assert 6*q<=96-4*len(T)
  U=np.zeros((len(T),96),dtype=np.int8)
  for i,t in enumerate(T):U[i,4*t:4*t+4]=1
  M=U[1:]-U[0] if equal else U;assert np.array_equal(M[:,[4*t for t in (T[1:] if equal else T)]],np.eye(len(M),dtype=np.int8)) and len(M)==96-d
  terms=dict(dense=dense,transversal24=trans,support24=total,support6=32*q**4,pure=192+4*len(T));n=sum(terms.values());rec=dict(dimension=d,section='equal' if equal else 'zero',rows=T,terms=terms,total=n,baseline=base[d]['total'],gain=n-base[d]['total'],physical_permutations=maps,pattern_histogram=dict(hist),family_points=dict(family_points),projection_size=len(image),valid_images=good);reports.append(rec);print({k:rec[k] for k in ['dimension','total','gain','family_points']},flush=True)
args.output.write_text(json.dumps(dict(status='PASS',NR24_pairs=pairs,NR24_minimum_distance=minimum,records=reports),indent=2)+'\n')
