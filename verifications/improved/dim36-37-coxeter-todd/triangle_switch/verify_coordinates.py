"""Independent exact coordinate audit, without importing TriangleModel/search code.

Rebuild the base from a fresh rational shell enumeration and the source outer
code. Proposal triangle IDs are decoded from the fresh lexicographic norm-eight shell
and cross-checked against the explicit90-switch coordinate record.
Every pair involving a newly introduced point is tested by int64 products.
Unchanged pairs inherit the separately replayed exact source certificate.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib, json, math, time,sys,tempfile,shutil,subprocess
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SOURCE=ROOT
PROPOSAL=HERE/'proposal.json'
GEOMETRY=HERE/'switches.json'
start=time.monotonic()
def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()
c=json.loads((SOURCE/'completion/certificate.json').read_text())
p=json.loads(PROPOSAL.read_text())
assert p['format']=='ct_triangle_switch_delta_v1'
assert p['source_sha256']==digest(SOURCE/'completion/certificate.json')
assert __debug__, 'Run without -O'
meta=json.loads((HERE/'certificate.json').read_text());assert meta['proposal_sha256']==digest(PROPOSAL)
with tempfile.TemporaryDirectory(prefix='ct_independent_source_') as td:
 t=Path(td);(t/'phase').mkdir();(t/'completion').mkdir();shutil.copy2(SOURCE/'phase/exact_lattice.py',t/'phase/exact_lattice.py')
 for name in ['verify.py','certificate.json']:shutil.copy2(SOURCE/'completion'/name,t/'completion'/name)
 subprocess.run([sys.executable,'-B',str(t/'completion/verify.py')],check=True,capture_output=True,text=True)
 source_check=json.loads((t/'completion/verification.json').read_text());assert source_check['all_checks_pass'] and source_check['certificate_sha256']==p['source_sha256']
G=np.array(c['gram_12'],np.int64); M=np.array(c['sublattice_matrix'],np.int64)
GG=np.kron(np.eye(3,dtype=np.int64),G)
assert np.array_equal(G,G.T) and np.array_equal(M@G@M.T,3*G)
phi=np.eye(12,dtype=np.int64)-M
assert np.array_equal(phi@G@phi.T,G) and np.array_equal(phi@phi@phi,np.eye(12,dtype=np.int64))
first,second=c['label_first'],c['label_second']
assert sorted(first+second)==list(range(12))
powers=3**np.arange(6,dtype=np.int64)
def labels(v): return ((v[:,first]+v[:,second])%3)@powers

# Fresh rational LDL shell enumerator. No stored quotient or shell is trusted.
low=[[F(i==j) for j in range(12)] for i in range(12)]
diag=[]
for j in range(12):
    d=F(int(G[j,j]))-sum(low[j][k]**2*diag[k] for k in range(j))
    assert d>0
    diag.append(d)
    for i in range(j+1,12):
        low[i][j]=(F(int(G[i,j]))-sum(low[i][k]*low[j][k]*diag[k] for k in range(j)))/d
vec=[0]*12; shell=[]
def visit(j,left):
    if j<0:
        norm=8-left
        if norm:
            assert norm.denominator==1
            shell.append((tuple(vec),int(norm)))
        return
    offset=sum((low[i][j]*vec[i] for i in range(j+1,12)),F(0))
    sq=left/diag[j];rad=math.isqrt(sq.numerator//sq.denominator)+2
    center=math.floor(-offset)
    for value in range(center-rad,center+rad+1):
        cost=diag[j]*(F(value)+offset)**2
        if cost<=left:
            vec[j]=value;visit(j-1,left-cost)
visit(11,F(8))
assert Counter(n for _,n in shell)=={4:756,6:4032,8:20412}
V=np.array([v for v,n in shell],np.int64); norms=np.array([n for v,n in shell],np.int64)
assert np.array_equal(np.sum(V*(V@G),axis=1),norms)
lab=labels(V); weights=np.full(729,100,np.int64);weights[0]=0
for a,n in zip(lab,norms):weights[a]=min(weights[a],n)
cosets=[np.zeros((1,12),np.int64) if i==0 else V[(lab==i)&(norms==weights[i])] for i in range(729)]
assert Counter((int(weights[i]),len(cosets[i])) for i in range(729))=={(0,1):1,(4,3):252,(6,18):224,(8,81):252}
roots=np.flatnonzero(weights==4)
X0=np.array(c['root_phase_zero_vectors'],np.int64)
X=np.stack((X0,X0@phi,X0@phi@phi),axis=1)
assert np.array_equal(labels(X0),roots)
for r,x in zip(roots,X):assert set(map(tuple,x))==set(map(tuple,cosets[r]))
X=np.concatenate((X,np.zeros((1,3,12),np.int64)))

# Rebuild all 403704 source-base points independently, including the axial
# points of the zero outer coset: norm-four vectors times M have norm twelve.
outer=np.array(c['base_outer_cosets'],np.int64)
assert outer.shape==(1429,3) and np.all(outer[0]==0)
assert len(np.unique(outer,axis=0))==1429
blocks=[]; axial=V[norms==4]@M
for k in range(3):
    a=np.zeros((len(axial),3,12),np.int64);a[:,k]=axial;blocks.append(a.reshape(-1,36))
for row in outer[1:]:
    assert sum(weights[row])==12 and 0 in row
    blocks.append(np.array([np.concatenate(v) for v in product(*(cosets[x] for x in row))],np.int64))
B=np.concatenate(blocks)
assert B.shape==(403704,36)

def rowkeys(a):return set(map(tuple,a))
def validate_points(name,a,norm):
    a=np.asarray(a,np.int64)
    assert a.ndim==2 and a.shape[1]==36
    weighted=a@GG
    bound=int(np.abs(a).max(initial=0))*int(np.abs(weighted).max(initial=0))*36
    assert bound<2**63
    assert np.all(np.sum(a*weighted,axis=1)==norm),name
    assert len(np.unique(a,axis=0))==len(a),name
    return {'points':len(a),'squared_norm':norm,'unique':True,'integer_product_abs_bound':bound}
pointchecks={'source_base':validate_points('source base',B,12)}
# Triangle ID decoding is reconstructed from the documented owner/Y order.
Y=np.array(sorted(map(tuple,V[norms==8])),np.int64)
assert rowkeys(Y)==rowkeys(V[norms==8]) and len(Y)==20412
ylab=labels(Y); owners=[row for row in outer if sorted(map(int,weights[row]))==[0,4,8]]
assert len(owners)==756
flips=p['flipped_triangles']
assert all(type(t) is int and 0<=t<61236 for t in flips) and len(set(flips))==len(flips)==90
old=[];new=[];decoded=[]
for t in flips:
    owner,offset=divmod(t,81);row=owners[owner]
    i=int(np.flatnonzero(weights[row]==4)[0]);j=int(np.flatnonzero(weights[row]==8)[0]);k=3-i-j
    yy=Y[ylab==row[j]];assert len(yy)==81
    a=int(np.flatnonzero(roots==row[i])[0]);y=yy[offset]
    u=np.zeros((3,3,12),np.int64);w=u.copy()
    u[:,i]=X[a];u[:,j]=y;w[:,k]=-X[a];w[:,j]=y
    old.extend(u.reshape(3,36));new.extend(w.reshape(3,36))
    decoded.append({'triangle':t,'owner':owner,'Y_offset':offset,'old_root_block':i,'norm8_block':j,'new_negative_root_block':k,'root_index':a,'y':y.tolist()})
assert decoded==json.loads(GEOMETRY.read_text()), 'explicit switch decoding'
old=np.array(old,np.int64);C=np.array(new,np.int64)
oldkeys=rowkeys(old);Bkeys=rowkeys(B);Ckeys=rowkeys(C)
assert len(oldkeys)==len(Ckeys)==270 and oldkeys<=Bkeys and not(Ckeys&Bkeys)
U=np.array([v for v in B if tuple(v) not in oldkeys],np.int64)
assert len(U)==403434
assert len(set(t//81 for t in flips))>0
pointchecks['changed_base']=validate_points('changed base',C,12)

def expand(orbits,lift=False):
    assert all(len(row)==5 and all(type(x) is int for x in row) for row in orbits)
    a=np.array(orbits,np.int64)
    assert len(np.unique(a,axis=0))==len(a)
    assert np.all(a>=0) and np.all(a[:,:3]<(253 if lift else 252)) and np.all(a[:,3:]<3)
    out=[]
    for r in a:
        if lift:
            assert list(r[:3]).count(252)==1
            angles=(0,int(r[3]),int(r[4]));nz=[j for j in range(3) if r[j]!=252]
            assert angles[nz[0]]==angles[list(r[:3]).index(252)]==0
        for t in range(3):out.append(np.concatenate((X[r[0],t],X[r[1],(t+r[3])%3],X[r[2],(t+r[4])%3])))
    return np.array(out,np.int64)
source_orbits=c['phase_orbits']; source_set=set(map(tuple,source_orbits))
removed=p['remove_phase_orbits'];added=p['add_phase_orbits']
assert len(source_set)==len(source_orbits)==32661
assert len(set(map(tuple,removed)))==len(removed)==1469
assert len(set(map(tuple,added)))==len(added)==1511
assert set(map(tuple,removed))<=source_set
assert not(set(map(tuple,added))&source_set)
kept=[v for v in source_orbits if tuple(v) not in set(map(tuple,removed))]
R=expand(kept);N=expand(added);L=expand(c['lift_phase_orbits'],True)
assert len(R)==93576 and len(N)==4533 and len(L)==1908
for name,a,norm in [('retained_equator',R,12),('new_equator',N,12),('lift_bodies',L,8)]:
    pointchecks[name]=validate_points(name,a,norm)
final36=np.concatenate((U,C,R,N))
pointchecks['full_36']=validate_points('full36',final36,12)
assert len(final36)==501813
assert c['points_36']==501687 and c['points_37']==505503
assert c['lift_heads']==[2,-2] and len(c['lift_phase_orbits'])==636
assert len(final36)+2*len(L)==505629
print('Fresh source construction, norms, uniqueness, and actual switches passed.',flush=True)

categories={}
def cross(name,a,b,limit=6,same=False):
    aw=a@GG;best=-10**9;count=0;witness=None
    bound=int(np.abs(aw).max(initial=0))*int(np.abs(b).max(initial=0))*36
    assert bound<2**63
    for s in range(0,len(b),2048):
        d=aw@b[s:s+2048].T
        if same:
            ii=np.arange(s,min(s+2048,len(a)));d[ii,ii-s]=-10**6
        loc=np.unravel_index(np.argmax(d),d.shape);value=int(d[loc])
        if value>best:best=value;witness=[int(loc[0]),int(s+loc[1])]
        assert value<=limit,(name,value,witness)
        count+=d.size
    if same:count-=len(a)
    categories[name]={'ordered_point_pairs':count,'max_dot':best,'limit':limit,'max_witness_indices':witness,'pass':True}
    print(name,categories[name],flush=True)

cross('NEW_EQUATOR_vs_UNCHANGED_BASE',N,U)
cross('CHANGED_BASE_vs_UNCHANGED_BASE',C,U)
cross('CHANGED_BASE_vs_CHANGED_BASE',C,C,same=True)
cross('CHANGED_BASE_vs_RETAINED_EQUATOR',C,R)
cross('CHANGED_BASE_vs_NEW_EQUATOR',C,N)
cross('NEW_EQUATOR_vs_RETAINED_EQUATOR',N,R)
cross('NEW_EQUATOR_vs_NEW_EQUATOR',N,N,same=True)
cross('CHANGED_BASE_vs_LIFT_BODIES',C,L)
cross('NEW_EQUATOR_vs_LIFT_BODIES',N,L)
cross('LIFT_BODIES_vs_LIFT_BODIES_DISTINCT',L,L,limit=2,same=True)
# The two heights are explicitly represented, so every duplicate and full
# latitude pair is checked, including opposite heads with identical bodies.
W=np.concatenate((np.c_[L,np.full(len(L),2,np.int64)],np.c_[L,np.full(len(L),-2,np.int64)]))
H=np.zeros((37,37),np.int64);H[:36,:36]=GG;H[36,36]=1
assert len(np.unique(W,axis=0))==3816
assert np.all(np.sum(W*(W@H),axis=1)==12)
d=(W@H)@W.T;np.fill_diagonal(d,-1000000);assert d.max()<=6
categories['ACTUAL_37D_LATITUDES_vs_LATITUDES']={'ordered_point_pairs':len(W)*(len(W)-1),'max_dot':int(d.max()),'limit':6,'pass':True}
full37=np.concatenate((np.c_[final36,np.zeros(len(final36),np.int64)],W))
assert len(np.unique(full37,axis=0))==len(full37)==505629
assert np.all(np.sum(full37*(full37@H),axis=1)==12)
assert source_check['all_checks_pass'] and source_check['certificate_sha256']==p['source_sha256']
result={'all_checks_pass':True,'points_36':len(final36),'points_37':len(full37),'source_points_36':c['points_36'],'source_points_37':c['points_37'],'gains':[len(final36)-c['points_36'],len(full37)-c['points_37']],
    'source_sha256':p['source_sha256'],'proposal_sha256':digest(PROPOSAL),'triangle_geometry_sha256':digest(GEOMETRY),'fresh_rational_shell_counts':dict(Counter(map(int,norms))),
    'actual_binary_triangle_flips':len(flips),'old_base_points_removed':len(old),'new_base_points_added':len(C),'old_phase_orbits_removed':len(removed),'new_phase_orbits_added':len(added),
    'retained_lift_body_orbits':len(c['lift_phase_orbits']),'lift_heights':[2,-2],'point_checks':pointchecks,'pair_categories':categories,
    'inherited_pairs':'All pairs wholly inside unchanged source subsets use fresh isolated completion/verify.py replay. Every pair involving a new base or equator point is explicitly checked above. All actual latitude pairs are explicitly rechecked.',
    'unchanged_source_proof_replayed':True,'source_base_reconstructed_from_fresh_shells':True,'integer_arithmetic_overflow_bounds_checked':True,'full_37_norm_and_uniqueness':True,
    'independence':'No imports from TriangleModel, gain verifier, search code, or source enumerator. Independent rational shell enumeration and complete source-base reconstruction. Norm8 vectors are freshly enumerated and lexicographically ordered; all90 decoded switches must equal the explicit coordinate record.',
    'elapsed_seconds':time.monotonic()-start}
(HERE/'actual_switches.json').write_text(json.dumps(decoded,indent=2)+'\n')
assert (len(final36),len(full37))==(meta['points_36'],meta['points_37'])
(HERE/'coordinate_verification.json').write_text(json.dumps(result,indent=2)+'\n')
np.savez_compressed(HERE/'verified_delta_coordinates.npz',old_base=old,new_base=C,old_equator=expand(removed),new_equator=N,lift_bodies=L,gram_36=GG)
print(json.dumps(result,indent=2),flush=True)
