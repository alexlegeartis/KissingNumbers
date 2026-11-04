"""Export the exact upper-bound proof with explicit clique and PSD witnesses."""
import json
from fractions import Fraction
import numpy as np
from model import DATA,P

z=np.load(DATA/'label_moments.npz');O=z['orbits'];nv=len(z['reps'])
old=json.loads((DATA/'label_farkas972.json').read_text())
cuts=json.loads((DATA/'outer_cuts.json').read_text())['cuts']
four=np.load(DATA/'four_cliques.npy');base=2*nv+len(z['cliques'])+len(four)
lookup={}
for d in range(1024):
    row=np.bincount(O[d,0],minlength=nv)
    for p in range(21):
        q=row.copy();q[O[d^P[p],0,p]]+=1
        lookup.setdefault(tuple(map(int,q)),(d,p))
terms=[]
for i,w in old['inequality_terms']:
    if i<nv:term=dict(kind='nonnegative',variable=i)
    elif i<2*nv:term=dict(kind='at_most_one',variable=i-nv)
    elif i<2*nv+len(z['cliques']):
        d,p=lookup[tuple(map(int,z['cliques'][i-2*nv]))]
        term=dict(kind='clique',word=d,line=p)
    elif i<base:raise AssertionError('unexpected four-clique proof term')
    else:
        cut=cuts[i-base]
        term=dict(kind='psd',character=int(z['characters'][cut['block']]),
                  vector=cut['vector'],divisor=cut['divisor'])
    term['multiplier']=w;terms.append(term)
eq=[];edgevars=np.flatnonzero(z['edges'])
for i,w in old['equality_terms']:
    if i==0:term=dict(kind='diagonal')
    elif i<=len(edgevars):term=dict(kind='edge',variable=int(edgevars[i-1]))
    else:term=dict(kind='cardinality',value=old['target'])
    term['multiplier']=w;eq.append(term)
report=dict(target=old['target'],upper_component=old['upper_component'],
            fixed_base=27720,components=4,family_total_upper=27720+4*old['upper_component'],
            inequalities=terms,equalities=eq,contradiction=old['contradiction'],
            scope='Full line-profile family; no periodicity or independent-hole restriction',
            geometry_source='https://hexagonmath.org/2610.00248v1')
(DATA/'upper971_certificate.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(dict(inequalities=len(terms),equalities=len(eq),upper_component=971)),flush=True)
