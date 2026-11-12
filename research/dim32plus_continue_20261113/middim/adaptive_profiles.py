"""Numerical discovery of adaptive P48 profile compression; not a proof.

Input populations are the actual, complete sqrt(3)E8 fibers documented in
the pinned P48 head package. Heads are expressed in units of 1/sqrt(3),
making many natural radii rational. Source vectors are scaled to norm 18.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import sys,json,subprocess
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).parent
sys.path.insert(0,str(ROOT/'research/dim32plus_20261113/middim'))
import affine_profile_search as old

def build(kind,k):
    V=old.normals(kind,k)
    G=3*V@V.T/4
    assert np.all((old.u@V.T)%4==0)
    p,inv=np.unique(old.u@V.T//4,axis=0,return_inverse=True)
    w=np.bincount(inv,weights=old.weights).astype(np.int64)
    z=p@np.linalg.cholesky(3*np.linalg.inv(G))
    q=np.sum(z*z,axis=1)
    assert int(w.sum())==52416000
    return dict(V=V,G=G,p=p,w=w,z=z,q=q)

def run(spec,den=24,qmax=6.,seconds=45,minweight=1):
    S=build(spec[0],int(spec[1:])); vertices=[]
    for i in range(len(S['p'])):
        q=S['q'][i]
        if q>qmax+1e-9 or q>=18-1e-9 or S['w'][i]<minweight: continue
        r=np.sqrt(q)
        heights=sorted(set([j/den for j in range(-int(r*den+1e-9),int(r*den+1e-9)+1)]+[-r,r]))
        for h in heights:
            a=np.sqrt((18-h*h)/(18-q))
            vertices.append((i,a,h,int(S['w'][i]),*S['z'][i]))
    stem=f'{spec}_d{den}_q{qmax:g}'+(f'_w{minweight}' if minweight>1 else '')
    source=HERE/f'{stem}.txt'
    with source.open('w') as f:
        f.write(f'{len(vertices)} {S["p"].shape[1]}\n')
        for row in vertices:f.write(' '.join(map(str,row))+'\n')
    print('MODEL',spec,'profiles',len(S['p']),'vertices',len(vertices),flush=True)
    subprocess.run([str(HERE/'head_search'),str(source),str(HERE/f'{stem}_best.json'),'100000',str(seconds)],check=True)
    r=json.loads((HERE/f'{stem}_best.json').read_text())
    r.update(kind=spec,dimension=49-int(spec[1:]),denominator=den,qmax_scaled=qmax,
      normals=S['V'].tolist(),gram=S['G'].tolist(),slices=[dict(profile=S['p'][vertices[j][0]].tolist(),
      head=vertices[j][2],scale=vertices[j][1],weight=vertices[j][3],q_scaled=float(S['q'][vertices[j][0]])) for j in r['vertices']],
      status='Numerical discovery only; geometry and population provenance require independent exact replay')
    (HERE/f'{stem}_result.json').write_text(json.dumps(r,indent=2))
    print(json.dumps(r,indent=2),flush=True)
    return r

if __name__=='__main__':
    run(sys.argv[1],int(sys.argv[2]) if len(sys.argv)>2 else 24,
        float(sys.argv[3]) if len(sys.argv)>3 else 6.,float(sys.argv[4]) if len(sys.argv)>4 else 45.,int(sys.argv[5]) if len(sys.argv)>5 else 1)
