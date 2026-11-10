"""Candidate-specific wrong claims and invalid sections must be rejected."""
from copy import deepcopy
import json
from pathlib import Path
from section_replay import validate_claim
ROOT=Path(__file__).resolve().parents[1]
assert __debug__
records=json.loads((ROOT/'data/sections.json').read_text())
receipts={r['dimension']:r for r in json.loads((ROOT/'data/section_receipt.json').read_text())}
cases={}
for original in records:
 actual=receipts[original['dimension']]['terms'];validate_claim(original,actual)
 for label in ['wrong_total','wrong_sparse','wrong_dense','wrong_dimension','duplicate_row','wrong_RS','wrong_pure','wrong_transversal','wrong_gain']:
  bad=deepcopy(original)
  if label=='duplicate_row':bad['rows'][1]=bad['rows'][0]
  elif label in ['wrong_total','wrong_dimension','wrong_gain']:bad[label[6:]]+=1
  else:bad['terms'][{'wrong_sparse':'support24','wrong_dense':'dense','wrong_RS':'support6','wrong_pure':'pure','wrong_transversal':'transversal24'}[label]]+=1
  try:validate_claim(bad,actual)
  except AssertionError:cases[f"{original['dimension']}_{label}"]='rejected'
  else:raise AssertionError(label)
# Explicitly reject blindly using the independent-eight-row dense formula in the octad section.
bad=deepcopy(records[4]);assert bad['dimension']==88
bad['terms']['dense']=6**8*8**4
try:validate_claim(bad,receipts[88]['terms'])
except AssertionError:cases['octad_wrong_uniform_formula']='rejected'
else:raise AssertionError('incorrect octad formula accepted')
print(json.dumps({'status':'PASS','positive_controls':7,'rejected_cases':cases}))
