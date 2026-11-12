"""Search body-dependent medium splits with exact final local gates.

For each tested medium body, remove its current point, keep every other
medium point, and construct the polytope of heads allowed by the large and
medium layers. Tiny deletions are recorded jointly. The numerical polytope
only proposes points: all retained interactions are then checked exactly.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json,sys,math,time
import numpy as np
from scipy.spatial import HalfspaceIntersection,QhullError
from medium43_pool import W,P,HERE,cutoff,exact_beta
sys.path.insert(0,str(W));import verify_families as vf

def scan(perfamily=5,seed=20261113):
    started=time.time();rng=np.random.default_rng(seed)
    g,_,_,_=vf.geometry();f=np.load(W/'configuration_43.npz');old=np.load(W/'explicit_points_43.npz')
    raw=np.load(W/'tiny_bodies.npz')['H'].astype(np.int64);prior=json.loads((P/'candidate.json').read_text())
    D=int(prior['head_denominator']);C=f['classes'];R=f['r3'];d=int(f['head_denominator'])
    Zi=f['head_numerators'].astype(object)*(D//d);Z=Zi.astype(float)/D;FN=R+np.sum(Z*Z,axis=1)
    FNi=[int(r)*D*D+sum(int(a)*int(a) for a in z) for r,z in zip(R,Zi)]
    purei=f['pure_directions'];pure=purei.astype(float)
    med=list(map(int,np.flatnonzero(R==12)));major=np.flatnonzero(R!=12)
    H=np.stack([np.load(P/'medium_families'/f'family_{int(C[j])}.npz')['H'].astype(np.int64) for j in med])
    mid={j:i for i,j in enumerate(med)};flat=H.reshape(-1,48);deleted=set(prior['delete'])
    valid=np.array([j*495+row not in deleted for j in med for row in range(495)]).reshape(len(med),495)
    added=prior['points'];AB=np.stack([H[mid[p['family']],p['row']] for p in added]);AZi=np.array([p['head_numerators'] for p in added],dtype=object);AZ=AZi.astype(float)/D
    AN=12+np.sum(AZ*AZ,axis=1);ANi=[p['normal'] for p in added]
    keep=np.array([i not in set(prior['delete_tiny']) for i in range(len(old['body_indices']))]);oldidx=np.flatnonzero(keep)
    TI=old['body_indices'][keep];TB=raw[TI];TZi=old['head_numerators'][keep].astype(object)*(D//int(old['head_denominator']))
    tmeta=[dict(kind='old_tiny',index=int(i)) for i in oldidx]
    for p in prior['new_tiny']:
        TB=np.r_[TB,[raw[p['body_index']]]];TZi=np.r_[TZi,[np.array(p['head_numerators'],dtype=object)]];tmeta.append(dict(kind='prior_new_tiny',point=p))
    TZ=TZi.astype(float)/D;TN=10+np.sum(TZ*TZ,axis=1);TNi=[10*D*D+sum(int(a)*int(a) for a in z) for z in TZi]
    chosen=[]
    for j in med:
        rows=np.flatnonzero(valid[mid[j]])
        take=rows if perfamily>=len(rows) else rng.choice(rows,size=perfamily,replace=False)
        chosen.extend((j,int(row)) for row in take)
    BB=np.stack([H[mid[j],row] for j,row in chosen]);bmed=exact_beta(BB,flat).reshape(len(chosen),len(med),495)
    btiny=exact_beta(BB,TB);badded=exact_beta(BB,AB)
    options=[];candidates_seen=0
    for pos,(j,row) in enumerate(chosen):
        B=BB[pos];b=bmed[pos].copy();b[~valid]=-32768;b[mid[j],row]=-32768;maxbeta=b.max(axis=1)
        dc=((g['digits'][C]-g['digits'][int(C[j])])%3)@g['powers'];ub=12+R-18+g['mu'][dc]//4
        vv=np.r_[ub[major],maxbeta,badded[pos]];zz=np.r_[Z[major],Z[med],AZ];nn=np.r_[FN[major],FN[med],AN]
        zzi=np.r_[Zi[major],Zi[med],AZi];nni=[FNi[k] for k in major]+[FNi[k] for k in med]+ANi
        planes=np.c_[vv/np.sqrt(12*nn),2*zz/np.sqrt(nn[:,None]),-np.ones(len(nn))]
        planes=np.r_[planes,np.c_[np.zeros(len(pure)),2*pure/np.linalg.norm(pure,axis=1)[:,None],-np.ones(len(pure))],[[1,0,0,0,-1],[-1,0,0,0,0]]]
        try:vertices=HalfspaceIntersection(planes,np.array([.01,0,0,0]),qhull_options='QJ1e-10').intersections
        except QhullError:continue
        heads=[tuple(map(int,Zi[j]))]
        for v in vertices:
            if np.linalg.norm(v)<=1+1e-8 or v[0]<=1e-7:continue
            h=np.sqrt(12)*v[1:]/v[0]
            if np.max(np.abs(h))*D>=2**62:continue
            z=tuple(map(int,np.rint(h*D)))
            if z not in heads:heads.append(z)
        good=[]
        for z in heads:
            h=np.array(z,dtype=float)/D;n=12*D*D+sum(a*a for a in z);nf=12+h@h
            bad=btiny[pos]>np.sqrt(nf*TN)-2*(TZ@h)+1e-8
            if sum(bad)>2:continue
            if any(int(beta)>cutoff(z,n,w,m,D) for beta,w,m in zip(vv,zzi,nni)):continue
            if any((v:=2*sum(a*int(b) for a,b in zip(z,p)))>0 and v*v>n*sum(int(b)*int(b) for b in p) for p in purei):continue
            cut=np.array([cutoff(z,n,w,m,D) for w,m in zip(TZi,TNi)],dtype=np.int64)
            blocked=np.flatnonzero(btiny[pos]>cut).tolist()
            if len(blocked)<=2:good.append(dict(head_numerators=list(z),normal=n,blockers=blocked))
        candidates_seen+=len(good)
        for a,p in enumerate(good):
            for q in good[:a]:
                blockers=sorted(set(p['blockers'])|set(q['blockers']))
                if len(blockers)>2:continue
                if cutoff(p['head_numerators'],p['normal'],q['head_numerators'],q['normal'],D)<24:continue
                option=dict(family=j,row=row,points=[p,q],tiny_blockers=blockers)
                options.append(option)
                print('OPTION',j,row,'tiny deletions',len(blockers),'count',len(options),flush=True)
        if pos%20==0:print('SCAN',pos,'/',len(chosen),'retained head candidates',candidates_seen,'options',len(options),'seconds',time.time()-started,flush=True)
        if options and not options[-1]['tiny_blockers']:
            print('POSITIVE SINGLE SPLIT',flush=True);break
    result=dict(status='Each option exactly locally verified; cross-option angles and inherited-source replay still required',head_denominator=D,
      options=options,tiny_metadata=tmeta,checked_bodies=chosen[:pos+1],candidate_heads=candidates_seen,seed=seed,seconds=time.time()-started)
    (HERE/f'split43_scan_{seed}_{perfamily}.json').write_text(json.dumps(result,indent=2));print('DONE',len(options),'seconds',result['seconds'],flush=True)
    return result
if __name__=='__main__':scan(int(sys.argv[1]) if len(sys.argv)>1 else 5,int(sys.argv[2]) if len(sys.argv)>2 else 20261113)
