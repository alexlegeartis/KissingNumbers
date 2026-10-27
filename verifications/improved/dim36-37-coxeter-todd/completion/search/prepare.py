"""Build exact optional refill pools in a separate scratch directory.

The search is heuristic; completion/verify.py and verify_coordinates.py certify
its results independently. No generated search files belong in the package.
"""
import argparse
import itertools
import json
import hashlib
from pathlib import Path
import sys
import numpy as np

PACKAGE=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(PACKAGE/'phase'))
from exact_lattice import enumerate_short


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('--mode',choices=('equator','lifts'),default='equator')
    parser.add_argument('--certificate',type=Path,default=PACKAGE/'completion/certificate.json')
    args=parser.parse_args()
    output=args.output.resolve()
    if output==PACKAGE or PACKAGE in output.parents:
        parser.error('use a scratch directory outside the verification package')
    output.mkdir(parents=True,exist_ok=True)
    cert=json.loads(args.certificate.read_text())
    G=np.array(cert['gram_12'],dtype=np.int64);M=np.array(cert['sublattice_matrix'],dtype=np.int64)
    powers=3**np.arange(6,dtype=np.int64)
    weights=np.full(729,100,dtype=np.int64);weights[0]=0
    for v,n in enumerate_short(G.tolist()):
        label=sum(((v[a]+v[b])%3)*p for a,b,p in zip(cert['label_first'],cert['label_second'],powers))
        weights[label]=min(weights[label],n)
    roots=np.flatnonzero(weights==4);assert len(roots)==252
    digits=np.arange(729)[:,None]//powers%3
    D=weights[((digits[:,None,:]-digits[None,:,:])%3)@powers]
    base=np.array(cert['base_outer_cosets'],dtype=np.int64)
    zero=np.array(cert['root_phase_zero_vectors'],dtype=np.int64)
    phi=np.eye(12,dtype=np.int64)-M
    T=np.zeros((253,253,3),dtype=np.int8)
    for t in range(3):T[:252,:252,t]=(zero@G)@(zero@np.linalg.matrix_power(phi,t)).T
    T.tofile(output/'T.bin')
    if args.mode=='equator':
        allowed=[]
        for i,j in ((0,1),(0,2),(1,2)):
            safe=np.ones((252,252),dtype=bool)
            for row in base[base[:,3-i-j]==0]:
                safe &= D[roots,row[i]][:,None]+D[roots,row[j]][None,:]>=8
            allowed.append(safe)
        cells=[]
        for a in range(252):
            for b in np.flatnonzero(allowed[0][a]):
                for c in np.flatnonzero(allowed[1][a]&allowed[2][b]):cells.append((a,int(b),int(c)))
        pool=np.array([cell+(p,q) for cell in cells
                       for p,q in itertools.product(range(3),repeat=2)],dtype=np.int32)
        selected=cert['phase_orbits']
    else:
        rows=[]
        for k in range(3):
            i,j=[t for t in range(3) if t!=k]
            safe=np.ones((252,252),dtype=bool)
            for row in base:
                safe &= D[0,row[k]]+D[roots,row[i]][:,None]+D[roots,row[j]][None,:]>=8
            forbidden=np.zeros((252,252,3),dtype=bool)
            for orbit in cert['phase_orbits']:
                phases=(0,orbit[3],orbit[4])
                forbidden[orbit[i],orbit[j],(phases[j]-phases[i])%3]=True
            for a,b in zip(*np.where(safe)):
                for p in range(3):
                    if forbidden[a,b,p]:continue
                    labels=[252]*3;labels[i],labels[j]=int(a),int(b)
                    phases=[0]*3;phases[j]=p
                    rows.append(labels+phases[1:])
        pool=np.array(rows,dtype=np.int32)
        if 'lift_phase_orbits' in cert:selected=cert['lift_phase_orbits']
        else:
            index={int(label):i for i,label in enumerate(roots)};selected=[]
            for cell in cert['lift_body_cosets']:
                nonzero=[i for i,x in enumerate(cell) if x]
                for p in range(3):
                    phases=[0]*3;phases[nonzero[1]]=p
                    selected.append([index[x] if x else 252 for x in cell]+phases[1:])
    pool.tofile(output/'pool.bin')
    np.array(selected,dtype=np.int32).tofile(output/'selected.bin')
    metadata=dict(mode=args.mode,pool_orbits=len(pool),selected_orbits=len(selected),
                  baseline_sha256=hashlib.sha256(args.certificate.read_bytes()).hexdigest(),
                  maximum_dot=6 if args.mode=='equator' else 2,
                  binary_format='native little-endian int32 rows of five; int8 T[253][253][3]')
    assert sys.byteorder=='little'
    (output/'pool.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(json.dumps(metadata))


if __name__=='__main__':main()
