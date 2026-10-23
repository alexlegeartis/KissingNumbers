#!/usr/bin/env python3
"""Exact replay of the 41D, 42D and 44D constructions. See README.txt."""
import argparse,contextlib,concurrent.futures,json,subprocess,time
from pathlib import Path
from collections import Counter
from fractions import Fraction as F
import numpy as np
from flint import fmpz_mat
import geometry
import verify_masses
ROOT=Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--jobs',type=int,default=4)
    ap.add_argument('--skip-recount',action='store_true',help='Inspect algebra using recorded counts; not a full population replay')
    args=ap.parse_args();start=time.monotonic();g=geometry.geometry()
    with (ROOT/'reports/mass_replay.log').open('w') as f,contextlib.redirect_stdout(f):verify_masses.main()
    print('Parent embedding and exact 11-design population identities: PASS',flush=True)
    ai=np.array((4*fmpz_mat(g['p']['A'].tolist()).inv()).tolist(),dtype=np.int64)
    assert np.array_equal(g['p']['A']@ai,4*np.eye(8,dtype=np.int64))
    assert not np.any(g['u']@ai%4)
    powers=3**np.arange(8,dtype=np.int64);ids=((g['u']@ai//4)%3)@powers
    oldnew=np.empty(6561,dtype=np.int64);oldnew[g['uclasses']]=ids
    assert len(set(oldnew))==6561
    assert np.array_equal(oldnew[g['uclasses']],ids)
    newold=np.argsort(oldnew);digits=(np.arange(6561)[:,None]//powers)%3
    mu=np.full(6561,1000,dtype=np.int64);np.minimum.at(mu,ids,g['unorms'])
    assert np.array_equal(mu[oldnew],g['mu12'])
    core={int(r['class_id']):int(r['count']) for r in json.loads((ROOT/'data/core_counts.json').read_text())}
    certs={dim:dict(np.load(ROOT/f'certificate_{dim}.npz')) for dim in (41,42,44)}
    needed=sorted({(int(c),int(r)) for m in certs.values() for c,r in zip(m['classes'],m['r3']) if r>=16 and not (r==16 and mu[c]==8)})
    # Regenerate a rigorously enclosing integer LDL ellipsoid for each core shell.
    GF=fmpz_mat(g['q'].tolist());GI=GF.inv();tr=sum((F(str(GI[i,i])) for i in range(40)),F(0))
    P=2**24;C=np.linalg.cholesky(g['q'].astype(float));L=C/np.diag(C)[None,:]
    M=np.rint(L*P).astype(np.int64);D=np.rint(np.diag(C)**2*P).astype(np.int64)
    assert np.all(D>0) and not np.any(np.triu(M,1)) and np.all(np.diag(M)==P)
    err=P**3*GF-fmpz_mat(M.tolist())*fmpz_mat(np.diag(D).tolist())*fmpz_mat(M.tolist()).transpose()
    eta=max(sum(abs(int(err[i,j])) for j in range(40)) for i in range(40))
    work=ROOT/'reports/count_inputs';work.mkdir(exist_ok=True)
    jobs=[]
    for c,r in needed:
        target=3*r;bound=target*(P**3+eta*tr);radius=(bound.numerator+bound.denominator-1)//bound.denominator
        center=(g['labels'][newold[c]]+1)%3-1
        path=work/f'{c}.input'
        with path.open('w') as f:
            f.write(f'40 {P} 3 1\n')
            for arr in [D[None,:],M,g['q']]:
                for row in arr:f.write(' '.join(map(str,row))+'\n')
            f.write(f'{c} {target} {radius}\n'+' '.join(map(str,center))+'\n')
        jobs.append((c,path))
    if not args.skip_recount:
        exe=ROOT/'reports/exact_count'
        subprocess.run(['g++','-O3','-std=c++17',str(ROOT/'exact_count.cpp'),'-o',str(exe)],check=True)
        def run(job):
            c,path=job;p=subprocess.run([str(exe),str(path)],capture_output=True,text=True,check=True)
            result=json.loads(p.stdout);assert result['class_id']==c and result['count']==core[c]
            (ROOT/'reports'/f'count_{c}.json').write_text(p.stdout)
            return c
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            for c in pool.map(run,jobs):print('Core population:',c,core[c],'PASS',flush=True)
    pool=np.load(ROOT/'data/dual40_min.npz')['H'].astype(np.int64)
    assert len(pool)==30240 and np.all(np.sum(pool*pool,axis=1)==2880)
    assert not np.any(pool@g['p']['B'].T)
    co=pool@g['b'].T@g['qi3'];assert not np.any(co%288);co//=288
    assert np.array_equal(co@g['b'],pool) and len(np.unique(pool,axis=0))==len(pool)
    tiny=Counter(int(oldnew[g['lookup'][v.astype(np.uint8).tobytes()]]) for v in co%3)
    assert len(tiny)==1920
    reports=[]
    for dim,m in certs.items():
        cl=m['classes'];r3=m['r3'];z=m['head_numerators'];den=int(m['head_denominator']);weights=m['weights']
        assert all(m[k].dtype.kind in 'iu' for k in ['classes','r3','head_numerators','weights'])
        assert den>0 and len(cl)==len(r3)==len(z)==len(weights)
        assert fmpz_mat(z.tolist()).rank()==dim-40
        for c,r,w in zip(cl,r3,weights):
            if r==16 and mu[c]==8:assert w==116235
            elif r>=16:assert int(w)==core[int(c)]
            elif r==14:assert mu[c]==16 and w==9180
            elif r==12:assert mu[c]==24 and w==495
            else:assert r==10 and mu[c]==32 and int(w)==tiny[int(c)]
        n=len(cl);N=[int(r)*den**2+sum(int(a)**2 for a in v) for r,v in zip(r3,z)]
        checked=equal=0
        for i in range(n):
            diff=(digits[cl[i]]-digits[cl[:i+1]])%3
            b=r3[i]+r3[:i+1]-18+mu[diff@powers]//4
            same=(cl[:i+1]==cl[i])&(r3[:i+1]==r3[i]);b[same]=2*r3[i];b[i]=2*r3[i]-18
            dots=z[:i+1]@z[i]
            for j in range(i+1):
                p=int(b[j])*den**2+2*int(dots[j])
                if p>0:
                    slack=N[i]*N[j]-p*p;assert slack>=0,(dim,i,j,slack);equal+=slack==0
                checked+=1
        rec=dict(dimension=dim,points=int(weights.sum()),families=n,pairs=checked,equalities=equal,
            status='ALGEBRA_PASS_RECORDED_COUNTS' if args.skip_recount else 'PASS',
            types={str(int(r)):dict(families=int(np.sum(r3==r)),points=int(weights[r3==r].sum())) for r in np.unique(r3)})
        reports.append(rec);print(json.dumps(rec),flush=True)
    (ROOT/'reports/results.json').write_text(json.dumps(reports,indent=2))
    print('COMPLETE',time.monotonic()-start,'seconds',flush=True)

if __name__=='__main__':main()
