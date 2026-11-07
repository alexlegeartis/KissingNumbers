#!/usr/bin/env python3
"""Reconstruct and verify the frozen 769852-point 39D configuration independently.
Python 3 standard library plus a C++17 compiler. No producer code is imported.
Run: python3 verify_independent.py [--compiler c++]
"""
import argparse, hashlib, itertools, json, os, pathlib, platform, subprocess, time
HERE = pathlib.Path(__file__).resolve().parent
EXPECTED = {
    'orbit_certificate.json': 'dee846d08a9516ca12a55fcf89bde154e529ba6302d7146487bed92d18a58300',
    'base_certificate.json': '01bed5e6c9a138e13358e5b1e4d2cdc565135fce1c79ee17b01e1a305c37bad1',
    'supports_3424.txt': '5e750c76c21aa9e350784624d367cc724a392eb8cb006826613d8dc82321d1cf',
}
def require(condition, message):
    if not condition:
        raise ValueError(message)
def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    start = time.monotonic()
    parser = argparse.ArgumentParser(); parser.add_argument('--compiler', default=os.environ.get('CXX', 'c++'))
    args = parser.parse_args()
    # Respect repository compute guidance: reserve the highest allowed CPU.
    if hasattr(os, 'sched_getaffinity'):
        cpus = sorted(os.sched_getaffinity(0))
        if len(cpus) > 1: os.sched_setaffinity(0, cpus[:-1])
    for name, sha in EXPECTED.items():
        require(digest(HERE/'inputs'/name) == sha, 'Input hash: '+name)
    cert = json.loads((HERE/'inputs/orbit_certificate.json').read_text())
    base = json.loads((HERE/'inputs/base_certificate.json').read_text())
    require(cert['dimension'] == base['dimension'] == 39, 'dimension')
    require(cert.get('partial_middle') == [], 'this audit requires the complete middle layer')
    L,U=cert['norm_lower'],cert['norm_upper']
    require(type(L) is int and type(U) is int and (L,U)==(99999996,100000004),'norm envelope')
    G=cert['group_words']; A=cert['representatives']; n=len(A)
    require(len(G)==32 and all(type(g) is int and 0<=g<2**39 for g in G),'group words')
    group=set(G);require(len(group)==32 and 0 in group,'group uniqueness and identity')
    require(all(g^h in group for g in G for h in G),'group XOR closure')
    require(n==10269,'representative count')
    require(all(len(a)==39 and all(type(x) is int and 0<abs(x)<=U for x in a) for a in A),'dense full support and integer coordinate bounds')
    norms=[sum(x*x for x in a) for a in A]
    require(all(L*L<=t<=U*U for t in norms),'norm bounds (Python integers)')
    require(39*U**2 < 2**63,'int64 dot safety')
    # Every nonidentity group mask flips at least one nonzero coordinate.
    require(all(g==0 or any((g>>k)&1 for k in range(39)) for g in G),'free action')
    raw=(HERE/'inputs/supports_3424.txt').read_text().splitlines()
    require(len(raw)==3424 and all(len(s)==39 and set(s)<=set('01') and s.count('1')==8 for s in raw),'source supports')
    supports=[tuple(k for k,c in enumerate(s) if c=='1') for s in raw]
    changed=[]
    for r in base['support_replacements']:
        index=r['remove_index'];require(supports[index]==tuple(r['remove']),'replacement old support')
        new=tuple(r['add']); require(len(new)==8 and len(set(new))==8 and tuple(sorted(new))==new and all(0<=k<39 for k in new),'replacement new support')
        supports[index]=new;changed.append(index)
    require(len(changed)==len(set(changed))==4,'replacement count')
    require(len(set(supports))==3424,'distinct supports')
    support_masks=[sum(1<<k for k in s) for s in supports]
    max_overlap=0; overlap_hist={}; paircount=0
    for i,s in enumerate(support_masks):
        for t in support_masks[i+1:]:
            overlap=(s&t).bit_count();require(overlap<=4,'middle-middle support overlap')
            max_overlap=max(max_overlap,overlap);paircount+=1
            overlap_hist[overlap]=overlap_hist.get(overlap,0)+1
    even=[s for s in range(256) if s.bit_count()%2==0]
    require(len(even)==128 and len(set(even))==128,'128 even patterns')
    max_same_support=max(32-8*(a^b).bit_count() for i,a in enumerate(even) for b in even[i+1:])
    require(max_same_support<=16,'middle-middle same support')
    omitted=set(); heavy=[]
    for star in base['heavy_stars']:
        c=star['center'];leaves=star['leaves']
        require(type(c) is int and 0<=c<39 and len(leaves)==4 and len(set(leaves))==4 and c not in leaves and all(type(k) is int and 0<=k<39 for k in leaves),'heavy star schema')
        omitted.update(tuple(sorted((c,b))) for b in leaves)
        for cs in (-1,1):
            for mask in range(16):
                if mask.bit_count()%2:continue
                v=[0]*39;v[c]=4*cs
                for bit,k in enumerate(leaves):v[k]=-2 if (mask>>bit)&1 else 2
                heavy.append(tuple(v))
    require(len(omitted)==14 and len(heavy)==len(set(heavy))==64,'star/deleted-pair count')
    roots=[]
    for a,b in itertools.combinations(range(39),2):
        if (a,b) in omitted:continue
        for u,v in itertools.product((-4,4),repeat=2):
            row=[0]*39;row[a]=u;row[b]=v;roots.append(tuple(row))
    require(len(roots)==len(set(roots))==2908,'roots count')
    require(not set(roots)&set(heavy),'sparse duplicate')
    require(all(sum(x*x for x in s)==32 for s in roots+heavy),'sparse norm 32')
    require(all(sum(bool(x) for x in s)==2 for s in roots) and all(sum(bool(x) for x in s)==5 for s in heavy),'sparse shape distinction')
    # Any root-middle dot is at most 2*4*2=16, including all group signs.
    # Middle vectors have exactly eight nonzero entries. Dense have 39.
    sparse=3424*128+len(roots)+len(heavy); dense=32*n
    require(sparse==441244 and sparse+dense==cert['claimed_size']==769852,'total count')
    build=HERE/'_build';build.mkdir(exist_ok=True)
    data=build/'audit_input.txt'
    with data.open('w') as f:
        print(n,len(G),len(supports),len(roots),len(heavy),L,U,file=f)
        print(*G,file=f)
        for rows in (A,supports,roots,heavy):
            for row in rows:print(*row,file=f)
    executable=build/'audit_cpp'
    cmd=[args.compiler,'-O3','-std=c++17','-Wall','-Wextra','-Wpedantic',str(HERE/'audit.cpp'),'-o',str(executable)]
    subprocess.run(cmd,check=True)
    compiler_version=subprocess.run([args.compiler,'--version'],capture_output=True,text=True,check=True).stdout.splitlines()[0]
    with (HERE/'direct_results.json').open('w') as out, (HERE/'direct_progress.log').open('w') as log:
        subprocess.run([str(executable),str(data)],stdout=out,stderr=log,check=True)
    direct=json.loads((HERE/'direct_results.json').read_text()); require(direct['passed'],'C++ result')
    report={
        'all_checks_passed':True,'dimension':39,'lower_bound':sparse+dense,
        'certificate_sha256':EXPECTED['orbit_certificate.json'],
        'input_sha256':EXPECTED,'independent_python_sha256':digest(pathlib.Path(__file__)),
        'independent_cpp_sha256':digest(HERE/'audit.cpp'),
        'dense_representatives':n,'group_order':len(G),'dense_points':dense,
        'middle_supports':len(supports),'middle_points':len(supports)*128,
        'retained_roots':len(roots),'heavy_points':len(heavy),'sparse_points':sparse,
        'support_pairs_checked':paircount,'support_intersection_histogram':overlap_hist,
        'max_support_intersection':max_overlap,'max_same_support_middle_dot':max_same_support,
        'deleted_root_pairs':sorted(omitted),'direct_integer_checks':direct,
        'distinctness':'Sparse shapes and complete even signs; dense orbit pairs have normalized dot <= 1/2, ruling out all repeated rays and stabilizers; full dense support separates dense from sparse.',
        'compiler':compiler_version,'compile_flags':['-O3','-std=c++17','-Wall','-Wextra','-Wpedantic'],
        'platform':platform.platform(),'elapsed_seconds':time.monotonic()-start,
        'maximum_absolute_coordinate':max(abs(x) for a in A for x in a),'safe_coordinate_bound':U,
        'scope':'Full finite certificate; no lattice, solver, binary-code, floating-point, or producer-verifier assumption. No optimality or worldwide-priority claim. No remote comparison was performed.'}
    (HERE/'independent_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
