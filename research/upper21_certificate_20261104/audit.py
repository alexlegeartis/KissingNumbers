"""Check the real 960-point control, aggregate the certificates, and hash artifacts."""
import hashlib,json,os
from fractions import Fraction as F
import numpy as np
from model import HERE,DATA,OLD

z=np.load(DATA/'label_moments.npz');O=z['orbits'];nv=len(z['reps'])
L=np.load(HERE.parent/'attack21_structured_20261104'/'data'/'recovery.npz')['labels'][0]
u=np.flatnonzero(L>=0);n=len(u);assert n==960
ids=O[u[:,None]^u,L[u,None],L[u][None,:]]
counts=np.bincount(ids.ravel(),minlength=nv);sizes=np.bincount(O.ravel(),minlength=nv)
f=[F(21*int(c),n*int(s)) for c,s in zip(counts,sizes)]
q=np.load(DATA/'outer_lp.npz')
def dot(row):return sum((int(a)*v for a,v in zip(row,f)),F(0))
for row,b in zip(q['A'],q['b']):assert dot(row)>=int(b)
for row,b in zip(q['E'],q['e']):assert dot(row)==int(b)
assert dot(q['row'])==960
control=dict(status='PASS',selected=n,inequalities_checked=len(q['A']),
             equalities_checked=len(q['E']),arithmetic='exact Fraction from actual ordered-pair counts')
(DATA/'baseline_control.json').write_text(json.dumps(control,indent=2)+'\n')
upper=json.loads((DATA/'upper_verification.json').read_text())
relax=json.loads((DATA/'relaxation_verification.json').read_text())
assert upper['status']==relax['status']=='PASS'
runs=[json.loads(p.read_text()) for p in sorted(DATA.glob('*_result.json'))]
summary=dict(dimension=21,scope=upper['scope'],component_lower=960,component_upper=971,
             family_lower=31560,family_upper=31604,new_lower_bound=False,
             optimality_at_960_proved=False,upper_certificate=upper,
             relaxation_limit=relax,baseline_control=control,numerical_runs=runs)
(DATA/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
paths=list(HERE.glob('*.py'))+list(HERE.glob('*.cpp'))+[HERE/'README.md',HERE/'.gitignore']
paths += [p for p in DATA.iterdir() if p.suffix in ('.json','.npz','.npy','.log') and p.name!='manifest.json']
paths += [OLD/'component21.npz',OLD/'baseline21.npz',HERE.parent/'attack21_structured_20261104'/'data'/'recovery.npz']
manifest={os.path.relpath(p,HERE).replace(os.sep,'/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
(DATA/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(dict(status='PASS',component_interval=[960,971],family_interval=[31560,31604],
                      optimality_at_960_proved=False,hashed_files=len(manifest))),flush=True)
