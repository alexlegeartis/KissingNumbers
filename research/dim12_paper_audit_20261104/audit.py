"""Independent integer audit of the paper's construction and search symmetry."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='1'
from collections import Counter
import hashlib,io,itertools,json,re,zipfile
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
FACTORS=[[(0,1),(2,3),(4,5)],[(0,2),(1,4),(3,5)],[(0,3),(1,5),(2,4)],
         [(0,4),(1,3),(2,5)],[(0,5),(1,2),(3,4)]]


def clebsch_coefficients():
    """Rows are (A+B sqrt(2))/12; no floating point decides validity."""
    signs=np.array([s for s in itertools.product([-1,1],repeat=5) if np.prod(s)==1],dtype=np.int64)
    aa=[];bb=[]
    for floor in [0,1,-1]:
        for s in signs:
            a=np.zeros(6,dtype=np.int64);b=np.zeros_like(a)
            if floor==0:a[5]=4*s[4];b[1:5]=4*s[:4]
            else:a[0]=6*floor;a[5]=-6*s[4];b[1:5]=-3*s[:4]
            aa.append(a);bb.append(b)
    for i in range(6):
        for sign in [1,-1]:
            a=np.zeros(6,dtype=np.int64);a[i]=12*sign;aa.append(a);bb.append(np.zeros_like(a))
    aa=np.array(aa);bb=np.array(bb);a=[];b=[]
    for block in range(2):
        for ai,bi in zip(aa,bb):
            av=np.zeros(12,dtype=np.int64);bv=np.zeros_like(av)
            av[6*block:6*block+6]=ai;bv[6*block:6*block+6]=bi;a.append(av);b.append(bv)
    for matching in FACTORS:
        for (i,j),(k,l) in itertools.product(matching,repeat=2):
            for s in itertools.product([-1,1],repeat=4):
                av=np.zeros(12,dtype=np.int64);av[[i,j,k+6,l+6]]=6*np.array(s)
                a.append(av);b.append(np.zeros_like(av))
    return np.array(a),np.array(b)


def leq(a,b,bound):
    d=a-bound
    return ((b==0)&(d<=0) | (b>0)&(d<0)&(2*b*b<=d*d)
            | (b<0)&((d<=0)|(d*d<=2*b*b)))


def riesz(x,s):
    norm=np.linalg.norm(x,axis=1,keepdims=True);z=x/norm;g=z@z.T
    d=np.maximum(2-2*g,1e-14);logits=-s/2*np.log(d);np.fill_diagonal(logits,-np.inf)
    peak=logits.max();w=np.exp(logits-peak);total=w.sum();w/=total
    value=(np.log(total/2)+peak)/s
    grad=2*(w/d)@z;grad-=np.sum(grad*z,axis=1,keepdims=True)*z
    return float(value),grad/norm


def inspect_archive():
    rows=[]
    with zipfile.ZipFile(ROOT/'inputs'/'saved_48_systems.zip') as archive:
        for name in archive.namelist():
            if not name.endswith('.npy'):continue
            x=np.load(io.BytesIO(archive.read(name)),allow_pickle=False).astype(np.float64)
            assert x.shape==(48,6) and np.isfinite(x).all()
            x/=np.linalg.norm(x,axis=1,keepdims=True)
            g=x@x.T;values=g[np.triu_indices(48,1)]
            spec=np.linalg.eigvalsh(x.T@x)
            run=re.search(r'r(\d+).*b(\d+)',name)
            rows.append(dict(name=name,max_cos=float(values.max()),max_coordinate=float(abs(x).max()),
                             frame_spectrum=spec.tolist(),run=list(run.groups()) if run else None,
                             unit_violation=max(0.,float(values.max())-.5,float(abs(x).max())-.5)))
    spectra=Counter(tuple(np.round(r['frame_spectrum'],3)) for r in rows)
    return dict(arrays=len(rows),trajectory_labels=len(set(tuple(r['run']) for r in rows)),
                satisfying_constraints_to_1e8=sum(r['unit_violation']<=1e-8 for r in rows),
                satisfying_constraints_to_1e5=sum(r['unit_violation']<=1e-5 for r in rows),
                smallest_violation=min(r['unit_violation'] for r in rows),
                largest_violation=max(r['unit_violation'] for r in rows),
                approximate_frame_spectrum_groups=len(spectra),
                common_spectra=[dict(spectrum=list(k),count=v) for k,v in spectra.most_common(8)],
                rows=rows)


def main():
    a,b=clebsch_coefficients();assert a.shape==(840,12)
    assert np.all(np.sum(a*a+2*b*b,axis=1)==144)
    assert np.all(np.sum(a*b,axis=1)==0)
    pp=a@a.T+2*b@b.T;qq=a@b.T+b@a.T
    valid=leq(pp,qq,72);np.fill_diagonal(valid,True);assert valid.all()
    x=(a+np.sqrt(2)*b)/12
    signs=np.array(list(itertools.product([-1,1],repeat=12)),dtype=np.int64)
    pa=signs@a.T;pb=signs@b.T
    assert abs(pa).max()<64 and abs(pb).max()<64
    codes=(pa+64)*129+(pb+64)
    ordered=np.sort(codes,axis=1)
    unique,inverse,counts=np.unique(ordered,axis=0,return_inverse=True,return_counts=True)
    assert len(unique)==3 and sorted(counts.tolist())==[1024,1024,2048]
    parities=np.c_[np.prod(signs[:,1:6],axis=1),np.prod(signs[:,7:12],axis=1)]
    classes={}
    for c in range(3):
        ids=np.flatnonzero(inverse==c)
        classes[str(c)]=dict(count=len(ids),parity_pairs=np.unique(parities[ids],axis=0).tolist(),
                            representative=signs[ids[0]].tolist())
    # Sign flips with even parity on each five-coordinate equator preserve base.
    lookup={tuple(np.r_[ai,bi]):i for i,(ai,bi) in enumerate(zip(a,b))}
    rng=np.random.default_rng(260618984)
    flips=rng.choice([-1,1],size=12);flips[5]=np.prod(flips[1:5]);flips[11]=np.prod(flips[7:11])
    perm=np.array([lookup[tuple(np.r_[ai*flips,bi*flips])] for ai,bi in zip(a,b)])
    assert len(set(perm.tolist()))==840
    inv=np.argsort(perm);s=signs[18];xx=np.vstack([x,s/np.sqrt(12)])
    yy=np.vstack([x,(s*flips)/np.sqrt(12)]);fullperm=np.r_[inv,840]
    assert np.max(abs(yy-xx[fullperm]*flips))<1e-15
    # A signed-permutation equivariance check for the actual Adam recurrence.
    m1=np.zeros_like(xx);m2=np.zeros_like(yy);v1=np.zeros_like(xx);v2=np.zeros_like(yy)
    maximum=0.;trajectory=[]
    for step in range(1,11):
        _,g1=riesz(xx,8.);_,g2=riesz(yy,8.)
        m1=.9*m1+.1*g1;m2=.9*m2+.1*g2;v1=.999*v1+.001*g1*g1;v2=.999*v2+.001*g2*g2
        xx-=.001*(m1/(1-.9**step))/(np.sqrt(v1/(1-.999**step))+1e-8)
        yy-=.001*(m2/(1-.9**step))/(np.sqrt(v2/(1-.999**step))+1e-8)
        maximum=max(maximum,float(np.max(abs(yy-xx[fullperm]*flips))))
        trajectory.append(dict(step=step,coordinate_error=float(np.max(abs(yy-xx[fullperm]*flips))),
                               gradient_error=float(np.max(abs(g2-g1[fullperm]*flips)))))
    assert trajectory[0]['gradient_error']<1e-15
    # Check the theorem's displayed local inverse independently with integers.
    gamma9=np.array([[9,1,1,1,3],[1,9,1,1,3],[1,1,9,1,3],[1,1,1,9,3],[3,3,3,3,9]])
    inverse8=np.array([[9,0,0,0,-3],[0,9,0,0,-3],[0,0,9,0,-3],[0,0,0,9,-3],[-3,-3,-3,-3,12]])
    assert np.array_equal(gamma9@inverse8,72*np.eye(5,dtype=int))
    assert inverse8.sum()==24
    sources={p.name:dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
             for p in (ROOT/'inputs').iterdir() if p.is_file()}
    archive=inspect_archive();(ROOT/'archive_audit.json').write_text(json.dumps(archive,indent=2)+'\n')
    result=dict(paper='https://arxiv.org/abs/2606.18984v2',exact_clebsch840_pair_test=True,
                tested_sign_vectors=len(signs),inequivalent_initial_841_gram_histograms=len(unique),
                classes=classes,adam_roundoff_divergence_after_10_steps=maximum,
                adam_roundoff_trace=trajectory,
                local_inverse_and_condition4_exact=True,archive={k:v for k,v in archive.items() if k!='rows'},
                source_files=sources)
    (ROOT/'audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
