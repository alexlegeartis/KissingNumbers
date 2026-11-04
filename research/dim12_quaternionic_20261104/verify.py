"""Independent exact kissing test of the rays through binary64 coordinates.

For positive dot products, normalized cosine <= 1/2 is equivalent to
4*dot**2 <= norm_i**2 * norm_j**2. All these operations use Python integers.
No tolerance, floating point square root, optimizer, or Gram reconstruction
participates in the decision. This is a certificate of the supplied coordinates,
not of any optimality or impossibility claim.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np


def verify(path):
    path=Path(path)
    start=time.monotonic()
    x=np.load(path,allow_pickle=False) if path.suffix=='.npy' else np.loadtxt(path)
    x=np.asarray(x,dtype=np.float64)
    if x.ndim!=2 or x.shape[0]<2 or x.shape[1]!=12 or not np.isfinite(x).all():
        raise ValueError('Require a finite N by 12 array with N >= 2.')
    ratios=[[float(t).as_integer_ratio() for t in row] for row in x]
    exponent=max(d.bit_length()-1 for row in ratios for _,d in row)
    z=[[n << (exponent-(d.bit_length()-1)) for n,d in row] for row in ratios]
    norms=[sum(t*t for t in row) for row in z]
    if min(norms)<=0:
        raise ValueError('Zero vector has no direction.')
    violations=0
    first=None
    for i,row in enumerate(z):
        for j in range(i):
            dot=sum(a*b for a,b in zip(row,z[j]))
            if dot>0 and 4*dot*dot>norms[i]*norms[j]:
                violations+=1
                if first is None:
                    first=[j,i]
    # An independently normalized extended precision diagnostic; never the verdict.
    y=x.astype(np.longdouble)
    y/=np.sqrt(np.sum(y*y,axis=1,keepdims=True))
    peak=np.longdouble(-1)
    for i in range(1,len(y)):
        peak=max(peak,np.max(y[:i]@y[i]))
    return dict(path=str(path.resolve()),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                n=len(x),dimension=12,interpretation='exact binary64 rays, mathematically normalized',
                passed=violations==0,violating_pairs=violations,first_violating_pair_0based=first,
                max_cos_extended_precision=str(peak),integer_scale_power_of_two=exponent,
                elapsed_seconds=time.monotonic()-start)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('paths',nargs='+',type=Path)
    p.add_argument('--out',type=Path)
    a=p.parse_args()
    report=[verify(path) for path in a.paths]
    output=json.dumps(report,indent=2)+'\n'
    print(output,end='')
    if a.out:
        a.out.write_text(output)
    raise SystemExit(0 if all(r['passed'] for r in report) else 1)


if __name__=='__main__':
    main()
