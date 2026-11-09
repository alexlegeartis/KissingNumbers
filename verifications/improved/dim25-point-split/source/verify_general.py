#!/usr/bin/env python3
"""Standalone exact checker for rationally rotated Leech-circle constructions.
Input: only the JSON certificate. No construction code or saved shell is read.
The shell is independently rebuilt from the binary Golay polynomial.
"""
from pathlib import Path
from itertools import combinations, product
from collections import Counter
from hashlib import sha256
import argparse,json,math,re,time
import numpy as np
INT64_MAX=(1<<63)-1
SCALE=10**8

def require(ok,message):
    if not bool(ok):raise AssertionError(message)

def integer(x):
    require(type(x)is int or(type(x)is str and re.fullmatch(r"-?[0-9]+",x)),"Expected exact decimal integer")
    return int(x)

def ivec(x,n):
    require(type(x)is list and len(x)==n,"Wrong vector length")
    return tuple(integer(a)for a in x)

def dot(a,b):return sum(int(x)*int(y)for x,y in zip(a,b))

def independent_shell():
    """Enumerate a Golay/Leech shell without importing the supplied generator.

    The binary polynomial g has exponents {0,2,4,5,6,10,11}. We form all
    degree-at-most-11 binary multiples by integer XOR and extend by parity.
    The checked code weight distribution, linear construction, and cardinality
    imply distance >=8 and self-orthogonality. See PROOF.md for the elementary
    three-family shell pair bound.
    """
    polynomial = sum(1 << k for k in (0, 2, 4, 5, 6, 10, 11))
    words = []
    for message in range(1 << 12):
        word = 0
        for bit in range(12):
            if message & (1 << bit):
                word ^= polynomial << bit
        word |= (word.bit_count() & 1) << 23
        words.append(word)
    require(len(set(words)) == 4096, "Golay dimension is not 12")
    weights = Counter(word.bit_count() for word in words)
    require(weights == {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1},
            "Wrong extended Golay weight distribution")
    basis = [words[1 << bit] for bit in range(12)]
    require(all((x & y).bit_count() % 2 == 0 for x in basis for y in basis),
            "Golay basis is not self-orthogonal")

    shell = []
    # Family I: two coordinates of absolute value 4.
    for i, j in combinations(range(24), 2):
        for a, b in product((-4, 4), repeat=2):
            row = [0] * 24
            row[i], row[j] = a, b
            shell.append(tuple(row))
    n1 = len(shell)
    # Family II: an octad of absolute value 2, with an even number of minuses.
    for word in words:
        if word.bit_count() != 8:
            continue
        positions = [i for i in range(24) if word & (1 << i)]
        for negative in range(256):
            if negative.bit_count() % 2:
                continue
            row = [0] * 24
            for j, i in enumerate(positions):
                row[i] = -2 if negative & (1 << j) else 2
            shell.append(tuple(row))
    n2 = len(shell) - n1
    # Family III: a codeword sign vector with one sign-reversed triple.
    for word in words:
        signs = [-1 if word & (1 << i) else 1 for i in range(24)]
        for j in range(24):
            row = signs.copy()
            row[j] *= -3
            shell.append(tuple(row))
    n3 = len(shell) - n1 - n2
    require((n1, n2, n3) == (1104, 97152, 98304), "Wrong shell families")
    require(len(shell) == len(set(shell)) == 196560, "Wrong shell cardinality")
    require(all(dot(row, row) == 32 for row in shell), "Wrong shell norm")
    return shell, {str(k): v for k, v in sorted(weights.items())}


def quadratic_sign(a,b):
    """Exact sign of a+b*sqrt(3)."""
    a,b=int(a),int(b)
    if b==0:return (a>0)-(a<0)
    if a==0:return (b>0)-(b<0)
    sa,sb=(1 if a>0 else -1),(1 if b>0 else -1)
    if sa==sb:return sa
    d=a*a-3*b*b
    require(d!=0,"Nontrivial rational representation of sqrt(3)")
    return sa*(1 if d>0 else -1)

def field_sign(v):
    """Exact sign in basis 1,sqrt(2),sqrt(3),sqrt(6)."""
    a,b,c,d=map(int,v);sx=quadratic_sign(a,c);sy=quadratic_sign(b,d)
    if sx==0:return sy
    if sy==0:return sx
    if sx==sy:return sx
    delta=quadratic_sign(a*a+3*c*c-2*b*b-6*d*d,2*a*c-4*b*d)
    require(delta!=0,"Nontrivial sqrt(2) in Q(sqrt(3))")
    return sx*delta

def field_dot(x,y):
    a,b,c,d=x;e,f,g,h=y
    return [dot(a,e)+2*dot(b,f)+3*dot(c,g)+6*dot(d,h),
            dot(a,f)+dot(b,e)+3*dot(c,h)+3*dot(d,g),
            dot(a,g)+dot(c,e)+2*dot(b,h)+2*dot(d,f),
            dot(a,h)+dot(d,e)+dot(b,g)+dot(c,f)]

def checked_mm(a,b):
    require(a.dtype==np.int64 and b.dtype==np.int64,"Matrix dtype is not int64")
    bound=a.shape[1]*max(abs(int(a.min())),abs(int(a.max())))*max(abs(int(b.min())),abs(int(b.max())))
    require(bound<=INT64_MAX,"Potential signed int64 accumulation overflow")
    return a@b,bound

def reconstruct(cert,shell_set):
    q=integer(cert['rotation_half_tangent_denominator']);require(q>0,"Invalid rotation denominator")
    ps=list(map(integer,cert['rotation_half_tangent_numerators']))
    nc=3*sum(len(c['owners'])for c in cert['classes']);require(len(ps)==nc,"Rotation count mismatch")
    out=[];dens=[];labels=[];circle_ids=[];centres=[];circle_classes=[];removed=[];i=0;circle_id=0
    for k,cl in enumerate(cert['classes']):
        v=np.array(ivec(cl['shift'],24),dtype=object);require(dot(v,v)==48,"Circle shift norm is not48")
        for row in cl['owners']:
            u=np.array(ivec(row,24),dtype=object);other=u+v
            require(tuple(u)in shell_set and tuple(other)in shell_set,"Circle endpoints not in shell")
            require(dot(u,v)==-24,"Circle endpoint product is incorrect")
            removed.extend([tuple(u),tuple(other)]);centres.append(2*u+v);circle_classes.append(k)
            for j in range(3):
                # cos(theta),sin(theta) are represented with common denominator6.
                co=np.array([6,0,0,0]if j==0 else[0,0,-2,0],dtype=object)
                si=np.array([0,0,0,0]if j==0 else[0,0,0,2 if j==1 else-2],dtype=object)
                if k==0 and j==0:
                    co=np.array([0,3,1,0],dtype=object);si=np.array([3,0,0,-1],dtype=object)
                p=ps[i];rc=q*q-p*p;rs=2*q*p;rd=q*q+p*p
                co,si=rc*co-rs*si,rs*co+rc*si
                plus=co.copy();plus[0]+=6*rd;minus=-co.copy();minus[0]+=6*rd
                require(field_sign(plus)>=0 and field_sign(minus)>=0,"Circle foot is outside endpoint segment")
                # centre=(2u+v)/2, radius vector=cos(theta)*v/2,
                # last coordinate=2*sqrt(3)*sin(theta).
                a=np.zeros((4,25),dtype=object);a[:,:24]=co[:,None]*v[None,:];a[0,:24]+=(2*u+v)*(6*rd)
                a[:,24]=[12*si[2],12*si[3],4*si[0],4*si[1]]
                out.append(a);dens.append(12*rd);labels.append(k);circle_ids.append(circle_id);i+=1
            circle_id+=1
    for x in cert['small']:
        a=np.zeros((4,25),dtype=object);a[0]=ivec(x['numerator'],25);d=integer(x['denominator']);require(d>0,"Invalid small denominator");out.append(a);dens.append(d)
    for name in ['local_owners','new_owners']:
        for row in cert.get(name,[]):
            u=ivec(row,24);require(u in shell_set,"Deleted local owner is not in shell");removed.append(u)
    return out,dens,set(removed),np.array(labels),np.array(circle_ids),np.array(centres,np.int64),np.array(circle_classes)

def enclosure(a,d):
    p=10**80;roots=[p]+[math.isqrt(k*p*p)for k in[2,3,6]];lo=[];hi=[]
    for j in range(25):
        low=high=0
        for k in range(4):
            v=int(a[k,j]);r=roots[k];s=r+(k!=0);low+=v*(r if v>=0 else s);high+=v*(s if v>=0 else r)
        lo.append(low*SCALE//(d*p));hi.append(-((-high*SCALE)//(d*p)))
    return np.array(lo,np.int64),np.array(hi,np.int64)

def verify(path,output):
    started=time.time();raw=path.read_bytes();cert=json.loads(raw);shell,weights=independent_shell();shell_set=set(shell)
    co,den,removed,labels,circle_ids,centres,circle_classes=reconstruct(cert,shell_set);n=len(labels);total=len(co)
    keep=np.array([u not in removed for u in shell]);bulk=np.array(shell,np.int64)[keep];size=len(bulk)+total
    require(size==integer(cert['points']),"Certificate point count is incorrect")
    for i,(a,d)in enumerate(zip(co,den)):
        require(field_dot(a,a)==[32*d*d,0,0,0],f'Point{i} norm is not exactly32')
    print('Independent shell, deleted endpoints, all norms and count passed.',flush=True)
    theorem_pairs=0
    for k in range(len(cert['classes'])):
        C=centres[circle_classes==k]
        if len(C)>1:
            G,_=checked_mm(C,C.T);np.fill_diagonal(G,-1);require(int(G.max())<=16,"Same-class circle-centre bound failed")
            theorem_pairs+=9*len(C)*(len(C)-1)//2
    intervals=[enclosure(a,d)for a,d in zip(co,den)];L=np.array([x[0]for x in intervals]);U=np.array([x[1]for x in intervals]);error=np.max(U-L,axis=1)
    require(int(error.max())<=2,"Coordinate enclosure too wide")
    l1=np.sum(abs(L),axis=1);matrix_bound=25*int(np.max(abs(L)))**2;correction_bound=2*int(error.max())*int(l1.max())+25*int(error.max())**2
    require(matrix_bound+correction_bound<=INT64_MAX,"Potential exceptional-pair interval overflow")
    exact=screened=0
    for start in range(0,total,128):
        end=min(total,start+128);base,_=checked_mm(L[start:end],L.T)
        bound=base+error[start:end,None]*l1[None,:]+l1[start:end,None]*error[None,:]+25*error[start:end,None]*error[None,:]
        for row,i in enumerate(range(start,end)):
            js=np.arange(i+1,total)
            if i<n:
                within=js<n;ok=np.ones(len(js),bool);j=js[within];ok[within]=~((labels[j]==labels[i])&(circle_ids[j]!=circle_ids[i]));js=js[ok]
            screened+=len(js)
            for j in js[bound[row,js]>16*SCALE*SCALE]:
                val=field_dot(co[i],co[int(j)]);val[0]-=16*den[i]*den[int(j)]
                require(field_sign(val)<=0,f'Exceptional conflict{i},{j}');exact+=1
    require(screened+theorem_pairs==total*(total-1)//2,"Exceptional-pair coverage mismatch")
    print('Every exceptional pair checked by rigorous intervals, exact signs, or same-class theorem.',flush=True)
    bulk_exact=0;bulk_l1=np.sum(abs(bulk),axis=1)
    for i in range(n,total):
        require(not np.any(co[i][1:]),"Small point is not rational")
        base,bound_mm=checked_mm(bulk,L[i,:24,None]);require(bound_mm+int(bulk_l1.max())*int(error[i])<=INT64_MAX,"Potential bulk interval overflow")
        bounds=base[:,0]+bulk_l1*error[i]
        for b in bulk[bounds>16*SCALE]:
            require(dot(b,co[i][0,:24])<=16*den[i],f'Bulk conflict at small point{i-n}');bulk_exact+=1
    print('Every retained-shell/small pair checked; circle feet lie on deleted endpoint segments.',flush=True)
    result=dict(verified=True,certificate=str(path.resolve()),certificate_sha256=sha256(raw).hexdigest(),checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),dimension=25,point_norm_squared=32,maximum_allowed_distinct_inner_product=16,packing_size=size,independent_shell_size=len(shell),removed_shell_points=len(removed),retained_shell_points=len(bulk),classes=len(cert['classes']),circles=n//3,circle_points=n,rational_small_points=total-n,exceptional_points=total,exceptional_pairs=total*(total-1)//2,same_class_circle_pairs=theorem_pairs,interval_screened_pairs=screened,exact_exceptional_fallbacks=exact,small_bulk_pairs=(total-n)*len(bulk),exact_bulk_fallbacks=bulk_exact,interval_scale=SCALE,golay_weight_distribution=weights,floating_point_geometry_used=False,construction_or_other_verifier_imported=False,comparison_number_supplied_by_requester=197584,difference_from_requester_comparison=size-197584,record_claim=False,seconds=time.time()-started)
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True);return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('certificate',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args();verify(a.certificate,a.output)
