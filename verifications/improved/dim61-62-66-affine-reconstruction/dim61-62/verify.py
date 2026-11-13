from pathlib import Path
import json,sys,hashlib
import numpy as np
HERE=Path(__file__).parent;ROOT=HERE.resolve().parents[3];P=ROOT/'verifications/improved/dim61-67-compatible-sign-families'
sys.path.insert(0,str(P));from subcode67 import code,expand,span
if not __debug__:raise RuntimeError('Run without -O; assertions are essential.')
c=json.load(open(HERE/'certificate.json'));fam=c['family'];removed=sum(1<<i for i in fam['removed_coordinates'])
assert fam==json.loads((P/'data/short61_62.json').read_text()), 'baseline family mismatch'
z=np.load(ROOT/'verifications/improved/dim61-68-95-96-code-layers/source/dim61_96/code64_supports.npz')
parent=[int(s) for s in z['words'] if not int(s)&removed];gen=[sum(int(x)<<i for i,x in enumerate(row)) for row in z['generator']]
fixed={r['support']:code(r,True) for k in ('old','new') for r in fam[k]}
modified={r['support']:code(r,True) for r in c['modified_parent_classes']};assert not(set(modified)&set(fixed));assert set(modified)<=set(parent)
added=[code(r) for r in c['new_classes']];assert len(set(a['s'] for a in added))==len(added);assert not set(a['s'] for a in added)&(set(parent)|set(fixed))
assert all(not a['s']&removed for a in added)
for a in [*added,*modified.values()]:
 assert a['shift'].bit_count()%2==1
 assert all(a['s']&g==0 or (a['s']&g).bit_count()%2==0 for g in gen)
 assert 65535 in span(a['checks'])
allold={s:fixed.get(s,modified.get(s,dict(s=s))) for s in parent}
oldnew=[fixed[r['support']] for r in fam['new']]
checks=0
for i,a in enumerate(added):
 for b in [*allold.values(),*oldnew,*added[:i]]:
  w=(a['s']&b['s']).bit_count();checks+=1
  if w<=8:continue
  assert w<=10 and 'dual' in b,(w,'unspecified')
  inter=a['s']&b['s'];common=[h for h in a['dual']&b['dual'] if h.bit_count()==8 and h&inter==h and ((a['phase']^b['phase'])&h).bit_count()%2]
  assert common,'no phase repair'
for a in modified.values():
 for b in oldnew:assert (a['s']&b['s']).bit_count()<=8
 # original supports unaffected irrespective signs, checked baseline
G=sum(a['count'] for a in added);assert G==c['gain']
report=dict(pass_=True,new_classes=len(added),gain=G,modified=len(modified),pair_checks=checks,totals={'61':88975146+G,'62':122539020+G},source_family_sha256=hashlib.sha256((P/'data/short61_62.json').read_bytes()).hexdigest(),scope='Exact new subcodes and affected pairs; inherits independently replayed existing family and transfer geometry')
(HERE/'verification.json').write_text(json.dumps(report,indent=2));print(report)
