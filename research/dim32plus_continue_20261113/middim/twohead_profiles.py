"""Two-head P48 profile discovery. Numerical proposals only."""
from adaptive_profiles import build,HERE
import numpy as np
import json,sys,subprocess

def run(spec,angles=24,qmax=4.5,seconds=90):
    S=build(spec[0],int(spec[1:]));vertices=[]
    fractions=[.25,.5,.75,1.]
    for i in range(len(S['p'])):
        q=S['q'][i]
        if q>qmax+1e-9 or q>=18-1e-9 or S['w'][i]<1:continue
        heads=[np.zeros(2)]
        if q>1e-9:
            for r in fractions:
                for j in range(angles):heads.append(r*np.sqrt(q)*np.array([np.cos(2*np.pi*j/angles),np.sin(2*np.pi*j/angles)]))
        for h in heads:
            a=np.sqrt((18-h@h)/(18-q))
            vertices.append((i,a,int(S['w'][i]),*h,*S['z'][i]))
    stem=f'{spec}_two{angles}_q{qmax:g}'
    source=HERE/f'{stem}.txt'
    with source.open('w') as f:
        f.write(f'{len(vertices)} {S["p"].shape[1]} 2\n')
        for row in vertices:f.write(' '.join(map(str,row))+'\n')
    print('MODEL',spec,'vertices',len(vertices),flush=True)
    subprocess.run([str(HERE/'multihead_search'),str(source),str(HERE/f'{stem}_best.json'),'200000',str(seconds)],check=True)
    r=json.loads((HERE/f'{stem}_best.json').read_text())
    r.update(kind=spec,dimension=50-int(spec[1:]),qmax_scaled=qmax,
      normals=S['V'].tolist(),gram=S['G'].tolist(),slices=[dict(profile=S['p'][vertices[j][0]].tolist(),
      head=list(vertices[j][3:5]),scale=vertices[j][1],weight=vertices[j][2],q_scaled=float(S['q'][vertices[j][0]])) for j in r['vertices']],
      status='Numerical discovery only; geometry and population provenance require independent exact replay')
    (HERE/f'{stem}_result.json').write_text(json.dumps(r,indent=2))
    print(json.dumps(r,indent=2),flush=True)
    return r
if __name__=='__main__':run(sys.argv[1],int(sys.argv[2]) if len(sys.argv)>2 else 24,float(sys.argv[3]) if len(sys.argv)>3 else 4.5,float(sys.argv[4]) if len(sys.argv)>4 else 90.)
