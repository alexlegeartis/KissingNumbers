import json,sys,numpy as np
from pathlib import Path
from copy import deepcopy
from guards import validate_source,count_pattern
ROOT=Path(__file__).resolve().parents[1];P=ROOT/'data';sys.path.insert(0,str(ROOT.parent/'source'));from section_replay import validate_claim
assert __debug__
data=json.load(open(P/'nr_partition.json'));words=json.load(open(P/'nr_pair4.json'));NR=set(data['NR_words']);parts=[{x^t for x in NR} for t in data['translate_representatives']]
validate_source(data,words,NR,parts,words)
cases={}
def rejected(name,fn):
 try:fn()
 except (AssertionError,IndexError,TypeError):cases[name]='rejected'
 else:raise AssertionError('Accepted '+name)
for kind in ['generator','partition','translation','span','NR_word']:
 bad=deepcopy(data)
 if kind=='generator':bad['octacode_generator'][0][0]^=1
 elif kind=='partition':bad['partition'][0][0]^=1
 elif kind=='translation':bad['translate_representatives'][0]^=1
 elif kind=='span':bad['span_basis'][0]=0
 else:bad['NR_words'][0]^=1
 rejected('wrong_'+kind,lambda:validate_source(bad,words,NR,parts,words))
bad=words.copy();bad[0]=bad[1];rejected('duplicate_NR24_word',lambda:validate_source(data,bad,NR,parts,words))
z=np.load(ROOT.parent/'data/supports.npz');linear=[0]
for row in z['sign_generator']:
 x=sum(int(v)<<j for j,v in enumerate(row));linear +=[w^x for w in linear]
for equal,file in [(False,'nonlinear_patterns.json'),(True,'equal_nonlinear_patterns.json')]:
 patterns=json.load(open(P/file));key=next(k for k,p in patterns.items() if p['family']=='NR' and len(sum(p['positions'],[]))>=2);pat=patterns[key];count_pattern(key,pat,equal,linear,words)
 for kind in ['family','translation','count','duplicate_position','group_size']:
  bad=deepcopy(pat)
  if kind=='family':bad['family']='unsupported'
  elif kind=='translation':bad['translation']=1<<24
  elif kind=='count':bad['count']+=1
  else:
   groups=[g for g in bad['positions'] if g]
   if kind=='duplicate_position':groups[0][1]=groups[0][0]
   else:groups[0].append(groups[0][0])
  rejected(f'{equal}_{kind}',lambda:count_pattern(key,bad,equal,linear,words))
for rec in json.load(open(P/'verified.json'))['records']:
 terms=rec['terms'];validate_claim(rec,terms)
 for key in ['dense','support24','transversal24','support6','pure']:
  bad=deepcopy(rec);bad['terms'][key]+=1;rejected(f"{rec['dimension']}_{key}",lambda:validate_claim(bad,terms))
 bad=deepcopy(rec);bad['rows'][1]=bad['rows'][0];rejected(f"{rec['dimension']}_repeated_row",lambda:validate_claim(bad,terms))
print(json.dumps({'status':'PASS','nonlinear_rejections':len(cases),'cases':cases}))
