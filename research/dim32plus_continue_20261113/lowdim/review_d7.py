"""Independent explicit-profile audit of nested D7/D4 population certificates.

No imports from extension_model or verify_extension. The finite domain is built
coordinate by coordinate; D4 orbits are graph components of explicit generators;
moments and averaged objectives are summed over all individual profiles.
"""
from collections import Counter, deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, product
from math import factorial, gcd, prod
from pathlib import Path
import json
import time

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]
GAMMA=HERE.parent/'gamma'
SOURCE=REPO/'verifications/improved/dim68-69-gamma72-cross-sections/weyl_averaged/embedding.json'
PIN='5f180d5d6c737010083e72210e16b92926ec2e79ed9a0891d0960a018b2d7842'


def profiles(k):
    # This enumeration does not start with coordinate-permutation representatives.
    result=[]
    def visit(prefix,choices,norm,largest,second):
        if len(prefix)==k:
            result.append(tuple(prefix));return
        for x in choices:
            aa=abs(x);l,s=sorted((largest,second,aa),reverse=True)[:2]
            if norm+x*x<=128 and l+s<=8:
                visit(prefix+[x],choices,norm+x*x,l,s)
    visit([],range(-8,9,2),0,0,0)
    visit([],range(-7,8,2),0,0,0)
    assert len(result)==len(set(result))
    return result


def d4_components(domain):
    remaining=set(domain);index={};components=[]
    while remaining:
        seed=min(remaining);queue=deque([seed]);members=[];remaining.remove(seed)
        while queue:
            h=queue.popleft();members.append(h)
            adjacent=[(-h[0],)+h[1:]]
            for j in range(3):
                v=list(h);v[j],v[j+1]=v[j+1],v[j];adjacent.append(tuple(v))
            a,b,c,d=h
            numerator=[a+b+c+d,a+b-c-d,a-b+c-d,a-b-c+d]
            assert all(x%2==0 for x in numerator)
            adjacent.append(tuple(x//2 for x in numerator))
            for v in adjacent:
                assert v in domain
                if v in remaining:remaining.remove(v);queue.append(v)
        ci=len(components)
        for h in members:index[h]=ci
        components.append(members)
    return index,components


def matmul(a,b):
    return [[sum(x*y for x,y in zip(row,col))for col in zip(*b)]for row in a]


def invdet(a):
    n=len(a);m=[[F(x)for x in row]+[F(i==j)for j in range(n)]for i,row in enumerate(a)]
    determinant=F(1)
    for col in range(n):
        pivot=next(i for i in range(col,n)if m[i][col])
        if pivot!=col:m[pivot],m[col]=m[col],m[pivot];determinant=-determinant
        v=m[col][col];determinant*=v;m[col]=[x/v for x in m[col]]
        for i in range(n):
            if i!=col:
                factor=m[i][col];m[i]=[x-factor*y for x,y in zip(m[i],m[col])]
    return [row[n:]for row in m],determinant


def embedding():
    raw=SOURCE.read_bytes();assert sha256(raw).hexdigest()==PIN
    lattice=json.loads(raw)['catalogue_gram']
    normals=json.loads((GAMMA/'d7_embedding.json').read_text())['normals']
    assert len(normals)==7 and all(len(row)==72 and all(type(x)is int for x in row)for row in normals)
    E=[[1,-1,0,0,0,0,0],[0,1,-1,0,0,0,0],[0,1,1,0,0,0,0],
       [-1,0,0,1,0,0,0],[0,0,0,-1,1,0,0],[0,0,0,0,-1,1,0],
       [0,0,0,0,0,-1,1]]
    actual=matmul(matmul(normals,lattice),list(zip(*normals)))
    reference=[[4*x for x in row]for row in matmul(E,list(zip(*E)))]
    assert actual==reference
    Ei,det=invdet(E);assert abs(det)==2
    # Every D7 root has integer coefficients in this embedded normal basis.
    roots=[]
    for i,j in combinations(range(7),2):
        for a,b in product((-1,1),repeat=2):
            r=[0]*7;r[i]=a;r[j]=b
            coeff=matmul([r],Ei)[0];assert all(x.denominator==1 for x in coeff)
            ambient=matmul([coeff],normals)[0]
            assert all(x.denominator==1 for x in ambient)
            assert matmul(matmul([ambient],lattice),[[x]for x in ambient])[0][0]==8
            roots.append(tuple(8*x for x in r))
    assert len(roots)==len(set(roots))==84
    dual=[[2*x for x in row]for row in Ei]
    assert all(x.denominator==1 for row in dual for x in row)
    assert all(len({int(x)%2 for x in col})==1 for col in zip(*dual))
    return roots


def geometry(c):
    ps=list(map(tuple,c['profiles']));heights=list(map(F,c['heights']))
    assert len(ps)==len(set(ps))==len(heights)
    dot=lambda a,b:sum(F(x*y,16)for x,y in zip(a,b))
    q=[dot(h,h)for h in ps]
    assert all(t<8 for t in q)
    scale=[(8-t*t)/(8-n)for t,n in zip(heights,q)]
    assert all(x>0 for x in scale)
    gaps=[]
    for i,j in combinations_with_replacement(range(len(ps)),2):
        d=4-dot(ps[i],ps[j]);r=4-heights[i]*heights[j]
        assert d>=0 and r>=0
        g=r*r-scale[i]*scale[j]*d*d;assert g>=0,(i,j,str(g));gaps.append(g)
    if c['dimension']==69:assert all(t*t<=2 for t in heights)
    else:assert all(t==0 for t in heights)
    return gaps


def main():
    if not __debug__:raise RuntimeError('Run without -O')
    start=time.monotonic();normal_roots=embedding()
    domain4=set(profiles(4));index,components=d4_components(domain4)
    domain7=profiles(7)
    canonical=lambda h:tuple(sorted(abs(x)for x in h))
    counts=Counter(map(canonical,domain7))
    representatives=sorted(counts,key=lambda h:(sum(x*x for x in h),h))
    repindex={h:i for i,h in enumerate(representatives)}
    reports=[]
    for dim in (68,69):
        path=GAMMA/f'd7_compact{dim}.json';c=json.loads(path.read_text())
        assert c['population']==6218175600 and c['k']==7 and c['dimension']==dim
        assert c['representatives']==[list(h)for h in representatives]
        assert c['orbit_sizes']==[counts[h]for h in representatives]
        selection=list(map(tuple,c['profiles']));assert len(set(selection))==len(selection)
        hits=Counter(index[h]for h in selection)
        marginal={h:F(hits[index[h]],len(components[index[h]]))for h in domain4}
        assert sum(marginal.values())==len(selection)
        tags=[tuple(p)for p in c['tags']]
        assert all(len(p)==7 and all(type(x)is int and x>=0 and x%2==0 for x in p)
                   and sum(p)<=10 for p in tags)
        active=[[(i,p)for i,p in enumerate(pp)if p]for pp in tags]
        moments=[[0]*len(representatives)for pp in tags]
        objectives=[F(0)]*len(representatives)
        for h in domain7:
            col=repindex[canonical(h)]
            objectives[col]+=marginal[h[:4]]
            for row,powers in zip(moments,active):
                row[col]+=prod(h[i]**p for i,p in powers)
        known=sum((marginal.get(h[:4],F(0))for h in normal_roots),F(0))
        assert known==F(c['known_objective'])==0
        assert objectives==list(map(F,c['objective']))
        rows=[];rhs=[]
        for row,pp in zip(moments,tags):
            degree=sum(pp)//2
            numerator=6218175600*128**degree
            numerator*=prod(factorial(p)//(2**(p//2)*factorial(p//2))for p in pp)
            b=F(numerator,prod(72+2*j for j in range(degree)))
            b-=sum(prod(x**p for x,p in zip(h,pp))for h in normal_roots)
            rr=[x*b.denominator for x in row];bb=b.numerator;g=gcd(*rr,bb)
            rows.append([x//g for x in rr]);rhs.append(bb//g)
        assert rows==c['rows'] and rhs==c['rhs']
        dual=list(map(F,c['dual_exact']));assert len(dual)==len(rows)
        slack=[obj-sum(y*r[j]for y,r in zip(dual,rows))for j,obj in enumerate(objectives)]
        assert all(s>=0 for s in slack) and slack==list(map(F,c['slack_exact']))
        mean=sum(y*b for y,b in zip(dual,rhs))+known
        assert mean==F(c['exact_value'])
        ceiling=-((-mean.numerator)//mean.denominator)
        assert ceiling==c['exact_bound']
        bound=ceiling+(2 if dim==69 else 0)
        if dim==69:assert bound==c['with_two_poles']
        gaps=geometry(c)
        reports.append({'dimension':dim,'exact_mean':str(mean),'bound':bound,
                        'moments_directly_summed':len(tags),'geometric_pair_types':len(gaps),
                        'minimum_squared_gap':str(min(gaps)),
                        'minimum_positive_squared_gap':str(min(g for g in gaps if g>0)),
                        'minimum_dual_slack':str(min(slack)),
                        'certificate_sha256':sha256(path.read_bytes()).hexdigest()})
        print(reports[-1],flush=True)
    result={'status':'PASS','complete_profile_count':len(domain7),
            'signed_permutation_orbits':len(representatives),'known_D7_roots':len(normal_roots),
            'D4_profile_count':len(domain4),'D4_orientation_orbits':len(components),
            'normals_sha256':sha256((GAMMA/'d7_embedding.json').read_bytes()).hexdigest(),
            'reports':reports,'elapsed_seconds':time.monotonic()-start,
            'independence':'No producer/model/checker imports. Direct profile-domain traversal; D4 generator BFS; all moment/objective terms summed explicitly.',
            'published_premises':'Pinned Gamma72 catalogue identity, integrality, minimum8, shell population6218175600, spherical11-design.'}
    (HERE/'d7_independent_review.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
