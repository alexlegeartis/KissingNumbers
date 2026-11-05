"""Joint signed-octad exchange with independently selectable fiber parities.

The input and output are compressed signed-support constructions. Each chosen
full support may independently select its even or odd parity, with deletions;
new antipodal lines are selected individually. The resulting full geometry is
replayed by the existing independent signed-support verifier.
"""
from copy import deepcopy
from itertools import combinations
import numpy as np
from scipy.optimize import milp, Bounds, LinearConstraint
from scipy.sparse import coo_matrix
from ..core import Artifact, RunContext, Verification, check_params, integer
from .layers import verify_artifact as verify_layers


def _bad(a,b):
    s,v=a;t,w=b;o=s&t
    return abs(o.bit_count()-2*(o&(v^w)).bit_count())>4


def _lines(mask,parity=None):
    positions=[j for j in range(mask.bit_length()) if mask>>j&1]
    if len(positions)!=8:raise ValueError('support must have weight 8')
    return [(mask,sum(1<<positions[j+1] for j in range(7) if v>>j&1))
            for v in range(128) if parity is None or (v.bit_count()&1)==parity]


def signed_exchange(inputs,params,context):
    if len(inputs)!=1 or inputs[0].kind!='implicit_spherical_code' or inputs[0].data.get('backend')!='signed_support_layers':
        raise ValueError('signed exchange needs a signed-support construction')
    check_params(params,('candidate_supports','candidate_lines','flexible_supports','seconds','min_new_lines','return_neutral'))
    base=inputs[0];report=verify_layers(base,context)
    if report.status!='exact':raise ValueError('input must independently verify exactly')
    cert=base.data['certificate'];n=cert['dimension'];F=cert['full_even_supports'];E=[tuple(x) for x in cert['exceptional_antipodal_lines']]
    flex=params.get('flexible_supports',[])
    if not isinstance(flex,list) or any(type(x) is not int or x not in F for x in flex) or len(set(flex))!=len(flex):raise ValueError('flexible supports must be distinct incumbent full supports')
    context.limit('max_orbit_points',len(flex)*128)
    supports=params.get('candidate_supports',[]);raw=params.get('candidate_lines',[])
    if not isinstance(supports,list) or not isinstance(raw,list):raise ValueError('candidate lists required')
    if any(type(x) is not int or x<0 or x>=1<<n or x.bit_count()!=8 for x in supports):raise ValueError('invalid candidate support')
    context.limit('max_codewords',128*len(supports)+len(raw))
    Q=[]
    for x in supports:Q.extend(_lines(x))
    for row in raw:
        if not isinstance(row,list) or len(row)!=2:raise ValueError('invalid candidate line')
        s,v=row
        if type(s) is not int or type(v) is not int or s<0 or s>=1<<n or s.bit_count()!=8 or v<0 or v&s!=v or v&(s&-s):raise ValueError('invalid canonical antipodal line')
        Q.append((s,v))
    Q=sorted(set(Q));context.limit('max_codewords',len(Q));fixed=[s for s in F if s not in set(flex)]
    # A fixed complete parity fiber conflicts with every new line whose
    # different support meets it in >=5 coordinates. Same-support candidates
    # are excluded too: modifying that fiber requires listing it as flexible.
    Q=[x for x in Q if not any((x[0]&s).bit_count()>=5 for s in fixed)]
    nc=len(Q);context.limit('max_pairs',nc*(nc-1)//2+nc*(len(E)+len(flex)*128))
    def wrap(result_cert,status):
        provenance={'classification':'signed_fiber_search','solver_status':str(status)[:500],
                    'candidate_lines':len(Q),'flexible_fibers':len(flex)}
        out=Artifact('implicit_spherical_code',{'backend':'signed_support_layers','certificate':result_cert,'provenance':provenance}, {},Verification())
        out.verification=verify_layers(out,context);out.facts=out.verification.details['facts']
        return out
    min_new=integer(params.get('min_new_lines',0),'min_new_lines',0)
    neutral=params.get('return_neutral',False)
    if type(neutral) is not bool:raise ValueError('return_neutral must be boolean')
    if min_new>len(Q):raise ValueError('min_new_lines exceeds retained candidates')
    seconds=params.get('seconds',5)
    if type(seconds) not in (int,float) or not np.isfinite(seconds) or seconds<=0:raise ValueError('invalid solver time')
    context.limit('solver_seconds',seconds)
    if not Q:return wrap(deepcopy(cert),'no candidates after fixed-fiber filtering')
    # A matching from every candidate support to a distinct conflicting old
    # fiber proves no positive exchange: each support's selected lines are
    # bounded by deletions in its matched fiber. Preserve neutral-search mode.
    if not neutral:
        support_list=sorted({s for s,v in Q});match={}
        def augment(s,seen):
            for f,t in enumerate(flex):
                if f in seen or not 5 <= (s&t).bit_count() <= 7:continue
                seen.add(f)
                if f not in match or augment(match[f],seen):match[f]=s;return True
            return False
        if len(support_list)<=len(flex) and all(augment(s,set()) for s in support_list):
            return wrap(deepcopy(cert),'exact structural bound: no positive gain in this finite pool')
    # Variables: new lines, deleted existing exceptional lines, parity bit for
    # each flexible fiber, then deletion bits for its 64 even and 64 odd lines.
    nf=len(flex);ne=len(E);p0=nc+ne;d0=p0+nf;nv=d0+nf*128
    fiber=[(_lines(s,0),_lines(s,1)) for s in flex]
    rows=[];cols=[];values=[];upper=[]
    def constraint(terms,bound):
        r=len(upper)
        for c,v in terms:rows.append(r);cols.append(c);values.append(v)
        upper.append(bound)
    for i,j in combinations(range(nc),2):
        if _bad(Q[i],Q[j]):constraint(((i,1),(j,1)),1)
    by_support={}
    for i,(support,_) in enumerate(Q):by_support.setdefault(support,[]).append(i)
    for ids in by_support.values():constraint(tuple((i,1) for i in ids),64)
    # On two distinct octads meeting in 5, 6 or 7 positions, all 256
    # antipodal sign lines admit a partition into 64 incompatible quadruples.
    for s,t in combinations(by_support,2):
        if 5 <= (s&t).bit_count() <= 7:
            constraint(tuple((i,1) for i in by_support[s]+by_support[t]),64)
    for i,q in enumerate(Q):
        for j,e in enumerate(E):
            if _bad(q,e):constraint(((i,1),(nc+j,-1)),0)
        for f,(even,odd) in enumerate(fiber):
            for j,e in enumerate(even):
                if _bad(q,e):constraint(((i,1),(p0+f,-1),(d0+128*f+j,-1)),0)
            for j,e in enumerate(odd):
                if _bad(q,e):constraint(((i,1),(p0+f,1),(d0+128*f+64+j,-1)),1)
    # For two different octads with overlap 5, 6 or 7, the 64 old
    # antipodal lines and the 128 possible new lines split into 64 triangles:
    # each old line and one new line of each parity are mutually incompatible.
    # Consequently selected new lines cannot outnumber deleted old lines in
    # that fiber. This is valid for either chosen parity of the old fiber.
    for f,s in enumerate(flex):
        for support in {q[0] for q in Q if 5 <= (q[0]&s).bit_count() <= 7}:
            terms=[(i,1) for i,q in enumerate(Q) if q[0]==support]
            terms.extend((d0+128*f+j,-1) for j in range(128))
            constraint(terms,0)
    # All octads containing a fixed 6-set jointly support at most 64
    # antipodal lines. Orient each line at a fixed common coordinate; there
    # are 32 sign patterns on the six common coordinates. For each pattern,
    # lines on different octads conflict, and at most two lines on one octad
    # can coexist (their two remaining signs must both differ).
    def sixsets(support):
        positions=[j for j in range(n) if support>>j&1]
        return [sum(1<<j for j in group) for group in combinations(positions,6)]
    sixgroups={}
    for support in set(by_support)|set(flex)|{e[0] for e in E}:
        for h in sixsets(support):sixgroups.setdefault(h,set()).add(support)
    for h,supports in sixgroups.items():
        ids=[i for s in supports for i in by_support.get(s,[])]
        if not ids:continue
        fs=[f for f,t in enumerate(flex) if t in supports]
        es=[j for j,e in enumerate(E) if e[0] in supports]
        if len(supports)<2 and not fs and not es:continue
        terms=[(i,1) for i in ids]
        terms.extend((d0+128*f+j,-1) for f in fs for j in range(128))
        terms.extend((nc+j,-1) for j in es)
        constraint(terms,64-64*len(fs)-len(es))
    # Any signed octad code supported in one 9-set contains at most
    # 64 antipodal lines. Completing each zero coordinate in both ways
    # gives two antipodal edges per line; these edges form an induced
    # matching of the 9-cube, whose size is at most 2^(9-2)=128.
    supportsets=list(by_support);nine=set()
    for s,t in combinations(supportsets,2):
        if (s&t).bit_count()==7:nine.add(s|t)
    for s in supportsets:
        for t in flex:
            if (s&t).bit_count()==7:nine.add(s|t)
    for u in nine:
        ids=[i for i,q in enumerate(Q) if q[0]&u==q[0]]
        fs=[f for f,t in enumerate(flex) if t&u==t]
        es=[j for j,e in enumerate(E) if e[0]&u==e[0]]
        terms=[(i,1) for i in ids]
        terms.extend((d0+128*f+j,-1) for f in fs for j in range(128))
        terms.extend((nc+j,-1) for j in es)
        constraint(terms,64-64*len(fs)-len(es))
    for f in range(nf):
        for j in range(64):
            constraint(((d0+128*f+j,1),(p0+f,1)),1) # delete even only when p=0
            constraint(((d0+128*f+64+j,1),(p0+f,-1)),0) # odd only when p=1
    if min_new:constraint(tuple((i,-1) for i in range(nc)),-min_new)
    M=coo_matrix((values,(rows,cols)),shape=(len(upper),nv)).tocsr();obj=np.zeros(nv);obj[:nc]=-1;obj[nc:p0]=1;obj[d0:]=1
    result=milp(obj,integrality=np.ones(nv),bounds=Bounds(0,1),constraints=LinearConstraint(M,-np.inf,upper),options={'time_limit':float(seconds)})
    if result.x is None or result.fun is None or (result.fun>=-.5 and not (neutral and result.fun<=.5)):
        return wrap(deepcopy(cert),result.message)
    sol=result.x
    if np.any(abs(sol-np.rint(sol))>1e-6):raise RuntimeError('fractional solver incumbent')
    pick=[Q[i] for i in range(nc) if sol[i]>.5];removed_e={j for j in range(ne) if sol[nc+j]>.5};newE=[list(e) for j,e in enumerate(E) if j not in removed_e]
    for f,s in enumerate(flex):
        parity=int(sol[p0+f]>.5);kept=fiber[f][parity]
        for j,line in enumerate(kept):
            if sol[d0+128*f+64*parity+j]<.5:newE.append(list(line))
    newE.extend(map(list,pick));gain=len(pick)-len(removed_e)-sum(int(sol[j]>.5) for j in range(d0,nv))
    if gain<0 or (gain==0 and not neutral):raise RuntimeError('solver gain mismatch')
    newcert=deepcopy(cert);newcert['full_even_supports']=[s for s in F if s not in set(flex)];newcert['exceptional_antipodal_lines']=newE
    newcert['expected_cardinality']=int(cert['expected_cardinality'])+2*gain
    return wrap(newcert,result.message)

OPERATORS={'search.signed_fiber_exchanges':signed_exchange}
INPUT_KINDS={'search.signed_fiber_exchanges':('implicit_spherical_code',)}
OUTPUT_KINDS={'search.signed_fiber_exchanges':('implicit_spherical_code',)}
