"""Self-contained exact certificate: 1,126,314 kissing points in dimension 40.

Requires Python, NumPy, and g++ with __int128. The default fast pair check
represents bounded integers exactly in IEEE binary32; --integer-pairs uses
int64 exclusively, and is substantially slower. Neither mode uses tolerances.
"""
import os
# Bound memory/thread use before importing NumPy. User settings take precedence.
os.environ.setdefault('OPENBLAS_NUM_THREADS','4')
os.environ.setdefault('MKL_NUM_THREADS','4')
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,subprocess,time
import numpy as np
ROOT=Path(__file__).resolve().parent

def require(ok,message):
    if not ok:raise RuntimeError(message)

def positive_definite(q):
    n=len(q);L=[[F(int(i==j)) for j in range(n)] for i in range(n)];D=[]
    for i in range(n):
        for j in range(i):
            L[i][j]=(F(int(q[i,j]))-sum((L[i][k]*D[k]*L[j][k]
                for k in range(j)),F(0)))/D[j]
        d=F(int(q[i,i]))-sum((L[i][k]**2*D[k] for k in range(i)),F(0))
        require(d>0,'Gram matrix is not positive definite');D.append(d)
    determinant=F(1)
    for d in D:determinant*=d
    return determinant

def check_minimum(q,compiler,work):
    e=np.load(ROOT/'enclosure.npz',allow_pickle=False)
    require(all(e[k].dtype.kind in 'iu' for k in e.files),'Noninteger enclosure')
    P=int(e['P']);D=e['D'].astype(object);M=e['M'].astype(object)
    J=e['inverse_times_three'].astype(object);Q=q.astype(object)
    require(0<P<2**40 and D.shape==(40,) and M.shape==(40,40),'Malformed enclosure')
    require(np.all((D>0)&(D<2**36)) and np.max(np.abs(M))<2**40,'Enclosure bounds')
    require(np.all(np.diag(M)==P) and not np.any(np.triu(M,1)),'Not triangular')
    require(J.shape==(40,40) and np.array_equal(Q@J,3*np.eye(40,dtype=object)),
            'Incorrect inverse')
    trace=F(sum(int(J[i,i]) for i in range(40)),3)
    E=P**3*Q-(M*D[None,:])@M.T
    eta=max(sum(abs(int(x)) for x in row) for row in E)
    bound=4*(P**3+eta*trace);radius=(bound.numerator+bound.denominator-1)//bound.denominator
    require(0<radius<2**120,'Radius out of bounds')
    path=work/'minimum.input'
    with path.open('w') as f:
        f.write(f'40 {P} {radius}\n')
        for a in [D[None,:],M,Q]:
            for row in a:f.write(' '.join(map(str,row))+'\n')
    exe=work/('exact_minimum.exe' if os.name=='nt' else 'exact_minimum')
    subprocess.run([compiler,'-O3','-std=c++17',str(ROOT/'exact_minimum.cpp'),'-o',str(exe)],check=True)
    run=subprocess.run([str(exe),str(path)],capture_output=True,text=True,check=True)
    result=json.loads(run.stdout)
    require(result['nonzero_vectors_of_norm_at_most_4']==0,'Lattice minimum below 6')
    return result

def max_product(a,b,q,same,integer):
    # q@b^T and the norm checks use int64. With the checked input bounds,
    # every absolute bilinear accumulation is < 40^2*128^2*2^20 < 2^45.
    qb=q@b.astype(np.int64).T
    accumulation_bound=40*int(np.max(np.abs(a)))*int(np.max(np.abs(qb)))
    if integer:
        aa=a.astype(np.int64);bb=qb
    else:
        require(np.finfo(np.float32).nmant==23,'IEEE binary32 required')
        require(accumulation_bound<2**24,'Products exceed exact binary32 range')
        aa=a.astype(np.float32);bb=qb.astype(np.float32)
    maximum=-10**9
    # Moderate tiles avoid allocating a gigabyte-sized product matrix.
    for j in range(0,len(b),2048):
        end=min(j+2048,len(b))
        for i in range(0,len(a),4096):
            stop=min(i+4096,len(a));dots=aa[i:stop]@bb[:,j:end]
            if same:
                lo=max(i,j);hi=min(stop,end)
                if lo<hi:
                    k=np.arange(lo,hi);dots[k-i,k-j]=-1000000
            maximum=max(maximum,int(dots.max()))
        if len(b)>10000:
            print(f'  checked {end}/{len(b)} columns; maximum {maximum}',flush=True)
    return maximum,accumulation_bound

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--compiler',default='g++')
    ap.add_argument('--integer-pairs',action='store_true')
    args=ap.parse_args();start=time.monotonic();work=ROOT/'replay';work.mkdir(exist_ok=True)
    q=np.array(json.loads((ROOT/'gram.json').read_text()),dtype=np.int64)
    require(q.shape==(40,40) and np.max(np.abs(q))<2**20,'Bad Gram matrix')
    require(np.array_equal(q,q.T) and not np.any(np.diag(q)%2),'Gram matrix not even integral')
    require(positive_definite(q)==6561,'Unexpected determinant')
    minimum=check_minimum(q,args.compiler,work)
    print('Exact rank-40 lattice minimum check: PASS.',flush=True)
    core=np.load(ROOT/'core.npz',allow_pickle=False)['coefficients']
    aux=np.load(ROOT/'auxiliary.npz',allow_pickle=False)
    groups=[core,aux['norm36'],aux['norm30']];sizes=[1084590,40590,1134];norms=[54,36,30]
    for a,size,norm in zip(groups,sizes,norms):
        require(a.shape==(size,40) and a.dtype==np.int8,'Bad coefficient array')
        require(len(np.unique(a,axis=0))==len(a),'Repeated vector')
        for i in range(0,len(a),8192):
            x=a[i:i+8192].astype(np.int64)
            require(np.all(np.sum((x@q)*x,axis=1)==norm),'Wrong vector norm')
    require(not np.any((core.astype(np.int16)-core[0])%3),'Core is not one affine coset')
    # For distinct core rows a,b, (a-b)/3 is a nonzero integer vector.
    # Thus (a-b)Q(a-b)^T >= 9*6=54, and aQb^T <= 27.
    print('Explicit counts, distinctness, norms and core coset: PASS.',flush=True)
    pairs=[]
    for i,j in [(1,1),(2,2),(1,2),(0,2),(0,1)]:
        print(f'Checking group pair {i},{j} ...',flush=True)
        value,bound=max_product(groups[i],groups[j],q,i==j,args.integer_pairs)
        require(value<=0 or 4*value*value<=norms[i]*norms[j],'Angle violation')
        # Cross-family raw norms differ; this strict bound excludes normalized overlap.
        if i!=j:require(value*value<norms[i]*norms[j],'Normalized families overlap')
        pairs.append(dict(groups=[i,j],maximum_integer_dot=value,
                          absolute_accumulation_bound=bound))
        print(f'Group pair {i},{j}: maximum {value}, PASS.',flush=True)
    result=dict(status='PASS',dimension=40,points=sum(sizes),
                previous_local_construction=1092072,improvement_over_previous=34242,
                published_benchmark=1092000,improvement_over_published=34314,
                group_sizes=sizes,coefficient_squared_norms=norms,
                core_pair_dot_bound=27,pair_checks=pairs,
                minimum_enumeration=minimum,
                pair_arithmetic='int64' if args.integer_pairs else 'exact bounded integers in IEEE binary32',
                seconds=round(time.monotonic()-start,3))
    files=['gram.json','enclosure.npz','core.npz','auxiliary.npz','exact_minimum.cpp','verify.py']
    result['sha256']={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files}
    (work/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
