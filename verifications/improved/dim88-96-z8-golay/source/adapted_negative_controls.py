import json
from copy import deepcopy
from pathlib import Path
from adapted_common import check_pattern,place_signs
from section_replay import validate_claim
ROOT=Path(__file__).resolve().parents[1]
assert __debug__
cases={}
def rejects(name,fn):
 try:fn()
 except (AssertionError,IndexError,TypeError):cases[name]='rejected'
 else:raise AssertionError('Accepted '+name)
for equal,file,key in [(False,'adapted_patterns.json','4,0'),(True,'adapted_equal_patterns.json','0,0,4,0,0')]:
 patterns=json.loads((ROOT/'data'/file).read_text());p=patterns[key];check_pattern(key,p,equal)
 for label in ['duplicate_position','outside_position','wrong_group_size']:
  bad=deepcopy(p)
  if label=='duplicate_position':bad['positions'][0][1]=bad['positions'][0][0]
  elif label=='outside_position':bad['positions'][0][0]=24
  else:bad['positions'][0].append(bad['positions'][1][0])
  rejects(f'{equal}_{label}',lambda:check_pattern(key,bad,equal))
 # Real physical cells: four selected pairs plus sixteen outside positions.
 cells=[0,1,4,5,8,9,12,13]+list(range(32,48));T=[0,1,2,3]
 place_signs(cells,T,p['positions'],equal)
 bad=deepcopy(p['positions']);bad[1]=bad[0]
 rejects(f'{equal}_nonbijective_physical_mapping',lambda:place_signs(cells,T,bad,equal))
 rejects(f'{equal}_missing_physical_cell',lambda:place_signs(cells[:-1],T,p['positions'],equal))
for equal,file in [(False,'adapted_zero.json'),(True,'adapted_equal.json')]:
 for rec in json.loads((ROOT/'data'/file).read_text()):
  rec=dict(rec,section='equal' if equal else 'zero');terms=rec['terms'];validate_claim(rec,terms)
  for key in ['dense','support24','transversal24','support6','pure']:
   bad=deepcopy(rec);bad['terms'][key]+=1
   rejects(f"{rec['dimension']}_{key}",lambda:validate_claim(bad,terms))
  bad=deepcopy(rec);bad['rows'][1]=bad['rows'][0]
  rejects(f"{rec['dimension']}_duplicate_row",lambda:validate_claim(bad,terms))
# Parity-incompatible equal rows cannot carry positive cardinality.
p={'count':1,'positions':[]}
rejects('mixed_parity_nonzero_claim',lambda:check_pattern('1,1,1,0,0',p,True))
print(json.dumps({'status':'PASS','adapted_rejections':len(cases),'cases':cases}))
