"""Independent coordinate audit, without TriangleModel or delta-helper imports.

Rebuilds changes solely from source certificate, actual original base coordinates,
and proposal. Exhaustively checks EVERY affected category of point pairs.
Unchanged/unchanged pairs rely on the separately replayed retained source proof.
"""
from pathlib import Path
import json,hashlib,time,numpy as np,sympy as sp,argparse,sys,subprocess,tempfile,shutil,itertools
P=Path(__file__).resolve().parent;SRC=P.parent;start=time.time()
parser=argparse.ArgumentParser();parser.add_argument('--write-coordinates',action='store_true');parser.add_argument('--prefix',default='replay');args=parser.parse_args();prefix=args.prefix;certpath=SRC/'completion/certificate.json';proposalpath=P/'proposal.json';c=json.loads(certpath.read_text());p=json.loads(proposalpath.read_text());sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest();assert p['source_sha256']==sha(certpath)
meta=json.loads((P/'certificate.json').read_text());assert meta['proposal_sha256']==sha(proposalpath) and meta['source_sha256']==sha(certpath)
assert __debug__, 'Run without -O'
with tempfile.TemporaryDirectory(prefix='ct_retained_source_') as td:
 t=Path(td);(t/'phase').mkdir();(t/'completion').mkdir();shutil.copy2(SRC/'phase/exact_lattice.py',t/'phase/exact_lattice.py')
 for name in ['verify.py','certificate.json']:shutil.copy2(SRC/'completion'/name,t/'completion'/name)
 subprocess.run([sys.executable,'-B',str(t/'completion/verify.py')],check=True,capture_output=True,text=True)
 source_report=json.loads((t/'completion/verification.json').read_text());assert source_report['all_checks_pass'] and source_report['certificate_sha256']==meta['source_sha256']
sys.path.insert(0,str(SRC/'phase'));from exact_lattice import enumerate_short
G=np.array(c['gram_12'],np.int64);GG=np.kron(np.eye(3,dtype=np.int64),G);assert np.array_equal(G,G.T);assert all(sp.Matrix(G[:k,:k].tolist()).det()>0 for k in range(1,13))
phi=np.eye(12,dtype=np.int64)-np.array(c['sublattice_matrix'],np.int64);assert np.array_equal(phi@G@phi.T,G) and np.array_equal(phi@phi@phi,np.eye(12,dtype=np.int64));R0=np.array(c['root_phase_zero_vectors'],np.int64);R=np.stack((R0,R0@phi,R0@phi@phi),axis=1);roots=R.reshape(-1,12);assert len(np.unique(roots,axis=0))==756;assert np.all(np.sum((roots@G)*roots,axis=1)==4)
first=np.array(c['label_first']);second=np.array(c['label_second']);pow3=3**np.arange(6)
def label(v):return ((v[:,first]+v[:,second])%3)@pow3
# Rebuild the entire source base from a fresh exact short-shell enumeration.
short=list(enumerate_short(G.tolist()));SV=np.array([v for v,n in short],np.int64);SN=np.array([n for v,n in short],np.int64);assert dict(zip(*np.unique(SN,return_counts=True)))=={4:756,6:4032,8:20412}
SL=label(SV);SW=np.full(729,999,np.int64);SW[0]=0;np.minimum.at(SW,SL,SN);cosets=[np.zeros((1,12),np.int64)if h==0 else SV[(SL==h)&(SN==SW[h])]for h in range(729)];blocks=[];M=np.array(c['sublattice_matrix'],np.int64);axial=SV[SN==4]@M
for j in range(3):
 b=np.zeros((756,3,12),np.int64);b[:,j]=axial;blocks.append(b)
for row in c['base_outer_cosets'][1:]:blocks.append(np.array(list(itertools.product(*(cosets[h]for h in row))),np.int64))
B=np.concatenate(blocks);assert B.shape==(403704,3,12)
# All block vectors and quotient minima are reconstructed from the fresh source base.
U,inv=np.unique(B.reshape(-1,12),axis=0,return_inverse=True);BI=inv.reshape(-1,3);norm=np.sum((U@G)*U,axis=1);assert set(norm)=={0,4,6,8,12};h=label(U);w=np.full(729,999,np.int64);np.minimum.at(w,h,norm);assert set(w)=={0,4,6,8}
rootlabels=label(R0);assert np.array_equal(rootlabels,np.flatnonzero(w==4));rootmap={int(a):i for i,a in enumerate(rootlabels)};Y=U[norm==8];Ylabels=label(Y);assert len(Y)==20412
owners=[]
for row in np.array(c['base_outer_cosets'],np.int64):
 if sorted(w[row])!=[0,4,8]:continue
 i=int(np.flatnonzero(w[row]==4)[0]);j=int(np.flatnonzero(w[row]==8)[0]);k=3-i-j;a=rootmap[int(row[i])];ys=np.flatnonzero(Ylabels==row[j]);assert len(ys)==81;owners.append((i,j,k,a,ys))
assert len(owners)==756
flips=p['flipped_triangles'];assert len(flips)==len(set(flips));old=[];new=[]
for t in flips:
 assert 0<=t<61236;i,j,k,a,ys=owners[t//81];y=Y[ys[t%81]];o=np.zeros((3,3,12),np.int64);n=o.copy();o[:,i]=R[a];o[:,j]=y;n[:,k]=-R[a];n[:,j]=y;old.extend(o.reshape(3,36));new.extend(n.reshape(3,36))
C=np.array(new,np.int64);O=np.array(old,np.int64);assert len(np.unique(C,axis=0))==len(np.unique(O,axis=0))==3*len(flips)
# Byte keys are exact integer-coordinate identity, not numerical proximity.
def rowbytes(a):return np.ascontiguousarray(a).view(np.dtype((np.void,a.shape[1]*a.dtype.itemsize))).ravel()
remove=np.isin(rowbytes(B.reshape(-1,36)),rowbytes(O));assert remove.sum()==3*len(flips);Bkeep=B[~remove].reshape(-1,36);Bfinal=np.vstack((Bkeep,C));assert len(np.unique(Bfinal,axis=0))==403704
# Look up changed blocks in the independently reconstructed coordinate dictionary.
umap={tuple(v):i for i,v in enumerate(U)};CI=np.array([[umap[tuple(v)]for v in row]for row in C.reshape(-1,3,12)],np.int64);BFI=np.vstack((BI[~remove],CI));assert np.array_equal(U[BFI].reshape(-1,36),Bfinal)
sourcephase=list(map(tuple,c['phase_orbits']));removed=list(map(tuple,p['remove_phase_orbits']));added=list(map(tuple,p['add_phase_orbits']));assert len(removed)==len(set(removed)) and len(added)==len(set(added));assert set(removed)<=set(sourcephase);finalphase=[v for v in sourcephase if v not in set(removed)]+added;assert len(finalphase)==len(set(finalphase))==len(sourcephase)-len(removed)+len(added)

def phase_indices(orbits):
 v=np.array(orbits,np.int64);angles=np.c_[np.zeros(len(v),np.int64),v[:,3:]];assert np.all((v[:,:3]>=0)&(v[:,:3]<252));return np.stack((3*v[:,:3]+angles,3*v[:,:3]+(angles+1)%3,3*v[:,:3]+(angles+2)%3),axis=1).reshape(-1,3)
EI=phase_indices(finalphase);NI=phase_indices(added);E=roots[EI].reshape(-1,36);N=roots[NI].reshape(-1,36);assert len(np.unique(E,axis=0))==3*len(finalphase);assert np.all(np.sum((E@GG)*E,axis=1)==12);assert np.all(np.sum((Bfinal@GG)*Bfinal,axis=1)==12)
report={'source_certificate_sha256':sha(certpath),'source_short_shell_counts':{str(n):int(np.sum(SN==n))for n in [4,6,8]},'unchanged_source_proof_replayed':True,'proposal_sha256':sha(proposalpath),'reconstruction_uses_triangle_api':False,'base_points':len(Bfinal),'equatorial_points':len(E),'changed_base_points':len(C),'new_equatorial_points':len(N),'removed_base_points':int(remove.sum()),'point_pair_checks':{}}
def record(name,count,maximum):report['point_pair_checks'][name]={'count':int(count),'maximum_inner_product':int(maximum),'passed':True};print(name,count,maximum,'seconds',time.time()-start,flush=True)
# 1. Every new equatorial point versus EVERY final base point. This category is
# intentionally explicit; no coset-admission or switch-requirement API is used.
RU=(roots@G@U.T).astype(np.int16);maximum=-999
for a in range(0,len(NI),16):
 ids=NI[a:a+16];dot=RU[ids[:,0,None],BFI[None,:,0]]+RU[ids[:,1,None],BFI[None,:,1]]+RU[ids[:,2,None],BFI[None,:,2]];maximum=max(maximum,int(dot.max()));assert dot.max()<=6
record('new_equator_vs_all_final_base',len(NI)*len(Bfinal),maximum)
# 2. Changed base versus every retained old base point.
maximum=-999;CG=C@GG
for a in range(0,len(Bkeep),8192):
 dot=CG@Bkeep[a:a+8192].T;maximum=max(maximum,int(dot.max()));assert dot.max()<=6
record('changed_base_vs_retained_base',len(C)*len(Bkeep),maximum)
dot=CG@C.T;assert np.all(np.diag(dot)==12);np.fill_diagonal(dot,-999);assert dot.max()<=6;record('changed_base_internal',len(C)*(len(C)-1)//2,dot.max())
# 3. Changed base versus every final equatorial point.
maximum=-999
for a in range(0,len(E),8192):
 dot=CG@E[a:a+8192].T;maximum=max(maximum,int(dot.max()));assert dot.max()<=6
record('changed_base_vs_all_final_equator',len(C)*len(E),maximum)
# 4. Every added equatorial point versus every final equatorial point.
RR=(roots@G@roots.T).astype(np.int16);emap={tuple(v):i for i,v in enumerate(EI)};selfcol=np.array([emap[tuple(v)]for v in NI]);maximum=-999
for a in range(0,len(NI),32):
 ids=NI[a:a+32];dot=RR[ids[:,0,None],EI[None,:,0]]+RR[ids[:,1,None],EI[None,:,1]]+RR[ids[:,2,None],EI[None,:,2]];dot[np.arange(len(ids)),selfcol[a:a+32]]=-999;maximum=max(maximum,int(dot.max()));assert dot.max()<=6
record('new_equator_vs_all_final_equator',len(NI)*len(E)-len(NI),maximum)
# Rank36 is certified by untouched axial source points and exact rank12 in each block.
for j in range(3):
 other=[k for k in range(3)if k!=j];b=Bfinal.reshape(-1,3,12);ax=b[np.all(b[:,other]==0,axis=(1,2)),j];assert sp.Matrix(ax.tolist()).rank()==12
report.update(points_36=len(Bfinal)+len(E),rank_36=36,all_norms_squared=12,threshold=6,all_final_points_distinct=True,source36_comparison=int(c['points_36']));assert report['points_36']==403704+3*len(finalphase)
# 37D: retain both height signs of every source lift body.
lift=np.array(c['lift_phase_orbits'],np.int64);RX=np.concatenate((R,np.zeros((1,3,12),np.int64)));la=np.c_[np.zeros(len(lift),np.int64),lift[:,3:]];L=np.stack([RX[lift[:,:3],(la+t)%3]for t in range(3)],axis=1).reshape(-1,36);assert len(L)==1908 and len(np.unique(L,axis=0))==1908;assert np.all(np.sum((L@GG)*L,axis=1)==8);assert sorted(c['lift_heads'])==[-2,2]
dot=CG@L.T;assert dot.max()<=6;record('changed_base_vs_lift_bodies',len(C)*len(L),dot.max());dot=N@GG@L.T;assert dot.max()<=6;record('new_equator_vs_lift_bodies',len(N)*len(L),dot.max());dot=L@GG@L.T;np.fill_diagonal(dot,-999);assert dot.max()<=2;record('distinct_lift_body_pairs',len(L)*(len(L)-1)//2,dot.max());report.update(points_37=report['points_36']+2*len(L),rank_37=37,lift_body_points=len(L),lift_heads=[-2,2],source37_comparison=int(c['points_37']));assert report['points_37']==report['points_36']+3816
report['gain_36']=report['points_36']-report['source36_comparison'];report['gain_37']=report['points_37']-report['source37_comparison'];report['all_affected_categories_checked_directly']=True;report['unchanged_pairs_rely_on_separately_replayed_retained_source_proof']=True;report['seconds']=time.time()-start
if args.write_coordinates:np.savez_compressed(P/(prefix+'_vectors36.npz'),coefficients=np.vstack((Bfinal,E)).astype(np.int16),gram=GG,squared_norm=12)
if args.write_coordinates:np.savez_compressed(P/(prefix+'_lifts37.npz'),body_coefficients=L.astype(np.int16),body_gram=GG,heads=np.array([-2,2],np.int16),squared_norm=12)
assert (report['points_36'],report['points_37'])==(meta['points_36'],meta['points_37'])
report['all_checks_pass']=True
(P/'verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True)
