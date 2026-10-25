"""Rationalize and exactly verify the collective 197584-point deformation.

Circle coordinates lie in Q(sqrt(2),sqrt(3)); other moved points are
rational points on the squared-radius-32 sphere. Integer outward bounds
screen products, with exact arithmetic for every inconclusive pair.
The standard Leech minimal-shell kissing property is the lattice input.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
import argparse,hashlib,json,math,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
from leech import build as leech


def qsign(a,b,d=3):
    a=int(a);b=int(b)
    if not b:return (a>0)-(a<0)
    if not a:return (b>0)-(b<0)
    sa=1 if a>0 else -1;sb=1 if b>0 else -1
    if sa==sb:return sa
    q=a*a-d*b*b;return sa*((q>0)-(q<0))


def sign(v):
    a,b,c,d=map(int,v);sx=qsign(a,c);sy=qsign(b,d)
    if not sx:return sy
    if not sy:return sx
    if sx==sy:return sx
    return sx*qsign(a*a+3*c*c-2*b*b-6*d*d,2*a*c-4*b*d)


def dot(x,y):
    a,b,c,d=x;e,f,g,h=y
    return [int(a@e+2*(b@f)+3*(c@g)+6*(d@h)),int(a@f+b@e+3*(c@h)+3*(d@g)),
            int(a@g+c@e+2*(b@h)+2*(d@f)),int(a@h+d@e+b@g+c@f)]


def reconstruct(cert):
    numerators=list(map(int,cert['rotation_half_tangent_numerators']));q=int(cert['rotation_half_tangent_denominator'])
    coefficients=[];denominators=[];removed=[];labels=[];centres=[];i=0
    for k,cls in enumerate(cert['classes']):
        v=np.array(cls['shift'],dtype=object);assert int(v@v)==48
        for row in cls['owners']:
            u=np.array(row,dtype=object);assert int(u@u)==32 and int(u@v)==-24
            removed.extend([tuple(map(int,u)),tuple(map(int,u+v))]);centres.append(np.array(2*u+v,dtype=np.int64))
            for j in range(3):
                co=np.array([6,0,0,0] if j==0 else [0,0,-2,0],dtype=object)
                si=np.array([0,0,0,0] if j==0 else [0,0,0,2 if j==1 else -2],dtype=object)
                if k==0 and j==0:co=np.array([0,3,1,0],dtype=object);si=np.array([3,0,0,-1],dtype=object)
                p=numerators[i];C=q*q-p*p;S=2*q*p;D=q*q+p*p
                co,si=co*C-si*S,si*C+co*S
                out=np.zeros((4,25),dtype=object);out[:,:24]=co[:,None]*v[None,:];out[0,:24]+=(2*u+v)*(6*D)
                out[:,24]=4*np.array([3*si[2],3*si[3],si[0],si[1]],dtype=object)
                coefficients.append(out);denominators.append(12*D);labels.append(k);i+=1
    for x in cert['small']:
        out=np.zeros((4,25),dtype=object);out[0]=list(map(int,x['numerator']));coefficients.append(out);denominators.append(int(x['denominator']))
    removed.extend(map(tuple,cert['local_owners']))
    removed.extend(map(tuple,cert.get('new_owners',[])))
    return coefficients,denominators,removed,np.array(labels),np.array(centres)


def enclosures(co,den,scale):
    precision=10**75;roots=[precision]+[math.isqrt(d*precision*precision) for d in [2,3,6]]
    lower=[];upper=[]
    for j in range(25):
        lo=hi=0
        for k in range(4):
            v=int(co[k,j]);r=roots[k];s=r+(k!=0)
            lo+=v*(r if v>=0 else s);hi+=v*(s if v>=0 else r)
        bottom=den*precision;lower.append(lo*scale//bottom);upper.append(-((-hi*scale)//bottom))
    return np.array(lower,np.int64),np.array(upper,np.int64)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--certificate',type=Path,default=ROOT/'collective_record25_certificate.json')
    ap.add_argument('--out',type=Path,default=ROOT/'collective_scaffold25_verification.json');args=ap.parse_args()
    started=time.time()
    cert=json.loads(args.certificate.read_text());co,den,removed,labels,centres=reconstruct(cert);n=len(labels)
    additions=len(cert.get('new_owners',[]))
    assert n==2913 and len(co)==3013+2*additions and len(set(removed))==1990+additions
    M=leech().astype(np.int64);lookup={tuple(u):i for i,u in enumerate(M)};keep=np.ones(len(M),bool)
    for u in removed:keep[lookup[u]]=False
    bulk=M[keep];assert len(bulk)+len(co)==cert['points']==197583+additions
    for c,d in zip(co,den):assert dot(c,c)==[32*d*d,0,0,0]
    # Different circles in one class clear for every phase, since each
    # centre product is at most 4 and each circular radius is sqrt(12).
    circle_groups=labels[::3]
    for k in range(3):
        C=centres[circle_groups==k];g=C@C.T;np.fill_diagonal(g,-100);assert g.max()<=16
    scale=10**8;intervals=[enclosures(c,d,scale) for c,d in zip(co,den)];L=np.array([x[0] for x in intervals]);U=np.array([x[1] for x in intervals]);error=np.max(U-L,axis=1);assert error.max()<=2
    l1=np.abs(L).sum(1);exact_pairs=0;screened=0
    for start in range(0,len(co),128):
        end=min(start+128,len(co));bound=L[start:end]@L.T+error[start:end,None]*l1[None,:]+l1[start:end,None]*error[None,:]+25*error[start:end,None]*error[None,:]
        for row,i in enumerate(range(start,end)):
            candidates=np.arange(i+1,len(co));needed=np.ones(len(candidates),bool)
            if i<n:
                within=candidates<n
                j=candidates[within];needed[within]=~((labels[j]==labels[i])&(j//3!=i//3))
            candidates=candidates[needed];screened+=len(candidates)
            for j in candidates[bound[row,candidates]>16*scale*scale]:
                v=dot(co[i],co[j]);v[0]-=16*den[i]*den[j]
                assert sign(v)<=0,('exceptional conflict',i,int(j))
                exact_pairs+=1
    print('exceptional products passed','exact fallbacks',exact_pairs,'seconds',round(time.time()-started,2),flush=True)
    # Circle feet lie on deleted endpoint segments and hence clear every
    # retained original minimal vector. Other moved points are rational.
    bulk_fallbacks=0;bl1=np.abs(bulk).sum(1)
    for j in range(n,len(co)):
        assert not np.any(co[j][1:])
        bounds=bulk@L[j,:24]+bl1*error[j]
        for b in bulk[bounds>16*scale]:
            value=sum(int(x)*int(y) for x,y in zip(b,co[j][0,:24]))
            assert value<=16*den[j],('bulk conflict',j,b.tolist())
            bulk_fallbacks+=1
    result=dict(status='EXACT VERIFICATION PASSED',points=cert['points'],bulk=len(bulk),circle_vertices=n,other_points=len(co)-n,
                exceptional_pairs_checked=screened,exact_exceptional_fallbacks=exact_pairs,exact_bulk_fallbacks=bulk_fallbacks,
                certificate_sha256=hashlib.sha256(args.certificate.read_bytes()).hexdigest(),seconds=time.time()-started,
                scope='Exact norm and product checks for the stated configuration. Standard Leech-shell kissing property assumed.')
    args.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':main()
