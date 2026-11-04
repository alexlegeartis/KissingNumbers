#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
from orbits17 import DATA

p=argparse.ArgumentParser()
p.add_argument('name')
a=p.parse_args()
if a.name == 'full17':
    z=np.load(DATA/'pool17.npz')
    n=len(z['vertices']); weights=np.ones(n,dtype=int); seed=z['seed'];fixed=4322;dim=17
else:
    z=np.load(DATA/(a.name+'.npz'))
    weights=z['weights'];n=len(weights);seed=z['seed'];fixed=int(z['fixed']);dim=int(z['dimension'])
with (DATA/(a.name+'.search')).open('w') as f:
    f.write(f'{n} {fixed} {dim}\n')
    f.write(' '.join(map(str,weights))+'\n')
    f.write(str(len(seed))+'\n'+' '.join(map(str,seed))+'\n')
print(a.name,n,len(seed),fixed+int(weights[seed].sum()))
