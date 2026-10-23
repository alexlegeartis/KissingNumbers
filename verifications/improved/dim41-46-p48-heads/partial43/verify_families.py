"""NumPy-only exact replay of the included kissing configurations.

Published inputs: P48n has minimum squared norm 6 and its 52,416,000 minimal
vectors form an 11-design. Large body populations use population_note.txt;
this script does not re-enumerate all their vectors.
"""
from pathlib import Path
import hashlib,json,time
import numpy as np
ROOT=Path(__file__).resolve().parent

def mm(a,b):
    a=np.asarray(a,dtype=np.int64);b=np.asarray(b,dtype=np.int64)
    assert a.shape[-1]*int(np.max(abs(a),initial=0))*int(np.max(abs(b),initial=0))<2**63
    return a@b

def inverse(a,multiple=1):
    v=np.rint(np.linalg.inv(a.astype(float))*multiple).astype(np.int64)
    assert np.array_equal(mm(a,v),multiple*np.eye(len(a),dtype=np.int64))
    return v

def geometry():
    p=np.load(ROOT/'parent_embedding.npz');d=np.load(ROOT/'body_complement.npz');g=np.load(ROOT/'class_geometry.npz')
    B=p['BB'];G=p['G'];b=d['b'];q=d['q'];Gi=inverse(G);Ui=inverse(p['U'])
    assert np.array_equal(mm(B,B.T),96*G)
    assert np.array_equal(mm(p['T'],B),p['B'])
    assert np.array_equal(mm(p['B'],p['B'].T),96*p['H'])
    assert np.array_equal(4*p['H'],3*mm(p['A'],p['A'].T))
    assert np.array_equal(mm(b,b.T),96*q) and not np.any(mm(b,p['B'].T))
    assert np.array_equal(mm(mm(p['U'],G),p['T'].T),np.r_[np.eye(8,dtype=np.int64),np.zeros((40,8),dtype=np.int64)])
    db=mm(b,B.T);assert not np.any(db%96)
    coeff=mm(mm(db//96,Gi),Ui);assert not np.any(coeff[:,:8]);inverse(coeff[:,8:])
    qi3=inverse(q,3)
    assert np.array_equal(mm(p['A'],g['ai']),4*np.eye(8,dtype=np.int64))
    u=g['u'];norm=np.sum(u*u,axis=1)
    assert np.all(u%2==u[:,:1]%2) and np.all(u.sum(1)%4==0)
    assert np.array_equal(norm,g['norm']) and np.all(norm<=32)
    ac=mm(u,g['ai']);assert not np.any(ac%4)
    ids=((ac//4)%3)@g['powers'];assert np.array_equal(ids,g['ids'])
    pair=mm(u,p['A'].T);assert not np.any(pair%4)
    lifts=mm(mm(pair//4,p['U'][:8]),B);dots=mm(lifts,b.T);assert not np.any(dots%96)
    co=mm(dots//96,qi3);normal3=3*lifts-mm(co,b)
    assert not np.any(mm(normal3,b.T))
    assert np.array_equal(np.sum(normal3*normal3,axis=1),72*norm)
    residues=np.empty((6561,40),dtype=np.int64)
    for i,c in enumerate(ids):residues[c]=co[i]%3
    assert np.array_equal(residues[ids],co%3) and len(np.unique(residues,axis=0))==6561
    assert np.array_equal((g['digits']@residues[g['powers']])%3,residues)
    mu=np.full(6561,9999,dtype=np.int64);np.minimum.at(mu,ids,norm)
    assert np.all(mu<9999) and np.array_equal(mu,g['mu']);mu//=4
    X=np.load(ROOT/'tiny_bodies.npz')['H'].astype(np.int64)
    assert len(X)==30240 and len(np.unique(X,axis=0))==30240
    assert np.all(np.sum(X*X,axis=1)==2880) and not np.any(mm(X,p['B'].T))
    tc=mm(mm(X,b.T),qi3);assert not np.any(tc%288);tc//=288
    assert np.array_equal(mm(tc,b),X)
    lookup={row.tobytes():i for i,row in enumerate(residues)}
    tiny=np.bincount([lookup[row.tobytes()] for row in tc%3],minlength=6561)
    assert np.array_equal(tiny,np.load(ROOT/'population_table.npz')['tiny'])
    return dict(g),residues,mu,tiny

def configuration(path,g,residues,mu,tiny):
    m=np.load(path);dim=int(m['dimension']);c=m['classes'];r=m['r3'];w=m['weights']
    z=m['head_numerators'].astype(object);D=int(m['head_denominator']);pure=m['pure_directions'].astype(object)
    assert D>0 and z.shape[0]==len(c)==len(r)==len(w)
    assert np.all(np.isin(r,[10,12,14,16,18])) and np.all(mu[c]%6==(18-r)%6)
    assert np.all(c[r==18]!=0) and np.all(mu[c[r==14]]==4) and np.all(mu[c[r==12]]==6)
    expected=np.where(r==18,1084590,np.where(r==16,116235-40*tiny[c],np.where(r==14,9180,np.where(r==12,495,tiny[c]))))
    assert np.array_equal(w,expected)
    assert np.array_equal(m['family_residues'],residues[c])
    A=m['span_projector_numerator'].astype(object);a=int(m['span_projector_denominator'])
    assert a>0 and np.array_equal(A,A.T) and np.array_equal(A@A,a*A) and np.trace(A)==(dim-40)*a
    assert np.array_equal(z@A,a*z) and np.array_equal(pure@A,a*pure)
    # A symmetric idempotent has rank equal to its trace. This is an exact
    # dimension certificate, even when the heads use eight ambient coordinates.
    N=r.astype(object)*D*D+np.sum(z*z,axis=1);checks=equalities=0
    for i in range(len(c)):
        diff=((g['digits'][c[i]]-g['digits'][c[:i+1]])%3)@g['powers']
        beta=(r[i]+r[:i+1]-18+mu[diff]).astype(object)
        beta[(c[:i+1]==c[i])&(r[:i+1]==r[i])]=2*int(r[i]);beta[-1]=2*int(r[i])-18
        doubled=beta*D*D+2*z[:i+1]@z[i];slack=N[i]*N[:i+1]-doubled*doubled
        assert np.all((doubled<=0)|(slack>=0)),('family angle',i)
        equalities+=int(np.sum((doubled>0)&(slack==0)));checks+=i+1
    pn=np.sum(pure*pure,axis=1);assert np.all(pn>0)
    for i,v in enumerate(pure):
        doubled=2*z@v;assert np.all((doubled<=0)|(doubled*doubled<=pn[i]*N)),('pure/family',i)
        doubled=2*pure[:i]@v;assert np.all((doubled<=0)|(doubled*doubled<=pn[i]*pn[:i])),('pure/pure',i)
        checks+=len(c)+i
    return dict(status='PASS',dimension=dim,points=int(w.sum())+len(pure),families=len(c),pure_points=len(pure),exact_angle_checks=checks,equalities=equalities,
                head_span_dimension=dim-40,layers={str(int(R)):dict(families=int(np.sum(r==R)),points=int(w[r==R].sum())) for R in np.unique(r)},
                configuration_sha256=hashlib.sha256(path.read_bytes()).hexdigest())

def main():
    start=time.monotonic();g,residues,mu,tiny=geometry();reports=[]
    for path in sorted(ROOT.glob('configuration_*.npz')):
        report=configuration(path,g,residues,mu,tiny);reports.append(report);print(json.dumps(report,indent=2),flush=True)
    assert reports
    result=dict(configurations=reports,geometry='Exact parent embedding, full rank-40 complement, gluing map, and 30240 tiny vectors replayed.',populations='Large shell populations inherited from population_note.txt; not freshly enumerated.',seconds=time.monotonic()-start)
    (ROOT/'replay.json').write_text(json.dumps(result,indent=2))

if __name__=='__main__':main()
