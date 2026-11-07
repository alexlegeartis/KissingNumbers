#!/usr/bin/env python3
"""Resolve global extrema into actual point pairs using Python integers."""
import itertools,json,pathlib,hashlib,os
P=pathlib.Path(__file__).resolve().parent
if hasattr(os, 'sched_getaffinity'):
    cpus=sorted(os.sched_getaffinity(0))
    if len(cpus)>1:os.sched_setaffinity(0,cpus[:-1])
raw=(P/'inputs/orbit_certificate.json').read_bytes()
assert hashlib.sha256(raw).hexdigest()=='dee846d08a9516ca12a55fcf89bde154e529ba6302d7146487bed92d18a58300'
c=json.loads(raw);b=json.loads((P/'inputs/base_certificate.json').read_text());r=json.loads((P/'direct_results.json').read_text())
A=c['representatives'];G=c['group_words']
def action(v,g):return [-x if (g>>k)&1 else x for k,x in enumerate(v)]
def pair(v,w):
    d=sum(x*y for x,y in zip(v,w));nv=sum(x*x for x in v);nw=sum(x*x for x in w)
    assert d<=0 or 4*d*d<=nv*nw
    return {'dot':d,'norm_squared_a':nv,'norm_squared_b':nw,'exact_squared_angular_slack':nv*nw-4*d*d}
i,j,g=r['dense_witness'];dense=pair(A[i],action(A[j],g));assert dense['dot']==r['max_dense_dot'];dense['labels']=[i,j,g]
supports=[tuple(k for k,x in enumerate(s) if x=='1') for s in (P/'inputs/supports_3424.txt').read_text().splitlines()]
for change in b['support_replacements']:
    assert supports[change['remove_index']]==tuple(change['remove']);supports[change['remove_index']]=tuple(change['add'])
i,s=r['middle_witness'];best=-10**30;middle_pair=None
for g in G:
    v=action(A[i],g)
    for mask in range(256):
        if mask.bit_count()%2:continue
        w=[0]*39
        for bit,k in enumerate(supports[s]):w[k]=-2 if (mask>>bit)&1 else 2
        d=sum(v[k]*w[k] for k in supports[s])
        if d>best:best=d;middle_pair=(v,w,g,mask)
assert best==2*r['max_dense_middle_half_dot']
v,w,g,mask=middle_pair;middle=pair(v,w);middle.update({'representative':i,'support_index':s,'group_mask':g,'even_sign_mask':mask,'sparse_vector':w})
omitted={tuple(sorted((star['center'],k))) for star in b['heavy_stars'] for k in star['leaves']}
rootpairs=[p for p in itertools.combinations(range(39),2) if p not in omitted]
best=-1;rw=None
for i,v in enumerate(A):
    av=list(map(abs,v))
    for x,y in rootpairs:
        d=4*(av[x]+av[y])
        if d>best:best=d;rw=(i,x,y)
assert best==r['max_dense_root_dot'];i,x,y=rw;w=[0]*39
for k in (x,y):w[k]=4 if A[i][k]>=0 else -4
root=pair(A[i],w);root.update({'representative':i,'sparse_vector':w})
# This shortcut only locates an actual witness; the global audit exhaustively
# checks all 21,030,912 dense/heavy pairs without this formula.
best=-1;hw=None
for i,v0 in enumerate(A):
    for g in G:
        v=action(v0,g)
        for si,star in enumerate(b['heavy_stars']):
            center=star['center'];leaves=star['leaves']
            signs=[-1 if v[k]<0 else 1 for k in leaves]
            if signs.count(-1)%2:signs[min(range(4),key=lambda t:abs(v[leaves[t]]))]*=-1
            d=4*abs(v[center])+2*sum(sign*v[k] for sign,k in zip(signs,leaves))
            if d>best:best=d;hw=(i,g,si,signs)
assert best==r['max_dense_heavy_dot'];i,g,si,signs=hw;star=b['heavy_stars'][si];v=action(A[i],g);w=[0]*39
w[star['center']]=4 if v[star['center']]>=0 else -4
for sign,k in zip(signs,star['leaves']):w[k]=2*sign
assert signs.count(-1)%2==0
heavy=pair(v,w);heavy.update({'representative':i,'group_mask':g,'star_index':si,'sparse_vector':w})
out={'passed':True,'arithmetic':'Python arbitrary-precision integers','dense_dense':dense,'dense_middle':middle,'dense_root':root,'dense_heavy':heavy,'scope':'Checks that every reported global extremum has an actual legal pair and a nonnegative exact squared angular slack; global maximality is proved by the separate complete C++ audit.'}
(P/'witness_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
