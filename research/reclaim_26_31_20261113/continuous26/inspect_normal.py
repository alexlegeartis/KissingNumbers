"""Rebuild normal-module inputs and inspect all fixed-exterior profiles."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--dimension',type=int,choices=[26,27],required=True)
    ap.add_argument('--out',type=Path,default=Path(__file__).parent);a=ap.parse_args();d=a.dimension
    root=Path(__file__).resolve().parents[3]
    src=root.parent/f'PackingStar/25D-31D/new_bounds_20261111/D{d}_{201180 if d==26 else 209862}_coordinates.npz'
    z=np.load(src);I=np.r_[[0,2,4,5,6,10,11,23],np.arange(24,d)]
    if d==26:
        assert int(z['rational_denominator'])==1
        X=z['rational_points'].astype(np.float64)
        X=np.vstack([X,z['quadratic_A']+np.sqrt(3)*z['quadratic_B']])
    else:X=z['points'].astype(np.float64)/int(z['coordinate_denominator'])
    J=[i for i in range(d) if i not in I];inside=np.all(X[:,J]==0,axis=1)
    Y=X[inside][:,I];B,w=np.unique(X[~inside][:,I],axis=0,return_counts=True)
    assert np.array_equal(B,np.rint(B)), 'Exterior profiles in this module must be integral'
    assert np.all(np.abs((Y*Y).sum(1)-32)<1e-10)
    a.out.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(a.out/f'normal{d}.npz',Y=Y,B=B,inside=np.flatnonzero(inside),indices=I)
    rows=[]
    for i,y in enumerate(Y):
        V=np.vstack([Y,B]);C=V[np.abs(V@y-16)<1e-9]
        rank=int(np.linalg.matrix_rank(C-np.outer(C@y,y)/32,tol=1e-8))
        rows.append([i,len(C),rank])
    norms=(B.astype(np.int64)**2).sum(1)
    profiles={}
    for n in np.unique(norms):
        values,counts=np.unique(w[norms==n],return_counts=True)
        profiles[str(int(n))]=dict(zip(map(str,values),map(int,counts)))
    report=dict(dimension=d,source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),source_count=len(X),sector=len(Y),
                profiles=len(B),profile_multiplicities_by_norm=profiles,contacts=rows,
                contact_scope='Numerical contact/rank diagnostics; profile coordinates and multiplicities are integer exact')
    (a.out/f'normal{d}_inspection.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='contacts'}))


if __name__=='__main__':main()
