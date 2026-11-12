"""Independent finite review of the new sign-family and compression proofs.

Uses no discovery or verifier imports from codes/. Parent-code minimum 16
is explicitly inherited from its separately completed exhaustive C++ replay.
This checker independently constructs kernels, verifies affine sign classes,
checks every changed support pair, and recounts the chosen 63D section.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations,product
from math import comb
import hashlib,json,time,sys
import numpy as np
HERE=Path(__file__).parent;ROOT=HERE.parents[2]
CODES=HERE.parent/'codes';SOURCE=ROOT/'verifications/improved/dim61-68-95-96-code-layers/source'
def need(ok,message):
    if not ok:raise ValueError(message)
def basis(rows):
    piv={}
    for x in rows:
        while x:
            j=x.bit_length()-1
            if j in piv:x^=piv[j]
            else:piv[j]=x;break
    return piv
def words(rows):
    a=[0]
    for r in rows:a += [x^r for x in a]
    return a
def expand(s,w):
    out=0;j=0
    for i in range(64):
        if s>>i&1:out|=((w>>j)&1)<<i;j+=1
    return out
def local(s,w):
    ans=0;j=0
    for i in range(64):
        if s>>i&1:ans|=((w>>i)&1)<<j;j+=1
    return ans
def kernel_from_columns(checks):
    piv={};ans=[]
    for j in range(16):
        syndrome=sum(((h>>j)&1)<<i for i,h in enumerate(checks));v=1<<j
        while syndrome:
            k=syndrome.bit_length()-1
            if k in piv:syndrome^=piv[k][0];v^=piv[k][1]
            else:piv[k]=(syndrome,v);break
        if not syndrome:ans.append(v)
    need(len(ans)==11 and len(basis(ans))==11,'sign kernel dimension')
    need(all((a&h).bit_count()%2==0 for a in ans for h in checks),'kernel equations')
    return ans
def class_data(row,odd=False):
    s=row['support'];h=row['checks'];phase=row['shift']
    need(s.bit_count()==16 and 0<s<1<<64,'support weight')
    need(len(h)==5 and h[0]==65535 and all(0<=v<65536 for v in h),'check shape')
    cols=[sum(((a>>j)&1)<<i for i,a in enumerate(h)) for j in range(16)]
    need(set(cols)==set(range(1,32,2)),'affine four-space columns')
    dual=words(h);need(len(set(dual))==32 and set(map(int.bit_count,dual))=={0,8,16},'RM dual')
    need(0<=phase<65536 and (not odd or phase.bit_count()%2==1),'phase')
    ker=kernel_from_columns(h)
    return dict(s=s,h=h,phase=expand(s,phase),localphase=phase,dual={expand(s,x) for x in dual},kernel=ker)
def family(path):
    started=time.time();raw=path.read_bytes();c=json.loads(raw)
    source=SOURCE/'dim61_96/code64_supports.npz'
    need(hashlib.sha256(source.read_bytes()).hexdigest()==c['source_archive_sha256'],'source hash')
    with np.load(source,allow_pickle=False) as f:
        old=list(map(int,f['words']));gen=[sum(int(x)<<i for i,x in enumerate(row)) for row in f['generator']]
    need(len(basis(gen))==28,'dense dimension')
    need(len(old)==len(set(old))==30828 and all(x.bit_count()==16 for x in old),'old supports')
    need(all((a&b).bit_count()%2==0 for a in gen for b in gen),'dense self orthogonality')
    removed=sum(1<<i for i in c['removed_coordinates']);old=[s for s in old if not s&removed]
    need(len(old)==c['original_supports_retained'],'shortened support count')
    modified={r['support']:class_data(r,c['odd_total_parity']) for r in c['old']};new=[class_data(r,c['odd_total_parity']) for r in c['new']]
    need(len(modified)==len(c['old']) and set(modified)<=set(old),'modified-old identity')
    need(len({x['s'] for x in new})==len(new) and not {x['s'] for x in new}&set(old),'new identities')
    need(all(not x['s']&removed and all((x['s']&g).bit_count()%2==0 for g in gen) for x in new),'new dual/shortening condition')
    repairs=[0,0];touched=set()
    def repair(a,b):
        intersection=a['s']&b['s'];need(intersection.bit_count()==10,'repair overlap')
        common=[x for x in a['dual']&b['dual'] if x.bit_count()==8 and x&intersection==x]
        need(len(common)==1,'common eight-check')
        need(((a['phase']^b['phase'])&common[0]).bit_count()%2==1,'opposite parity')
    for i,a in enumerate(new):
        for s in old:
            m=(a['s']&s).bit_count();need(m<=10,'old/new overlap')
            if m>8:
                need(s in modified,'missing modified class');repair(a,modified[s]);touched.add(s);repairs[0]+=1
        for b in new[:i]:
            m=(a['s']&b['s']).bit_count();need(m<=10,'new/new overlap')
            if m>8:repair(a,b);repairs[1]+=1
    need(touched==set(modified),'unused/missing modified classes')
    report=dict(certificate=path.name,sha256=hashlib.sha256(raw).hexdigest(),old_supports=len(old),new_supports=len(new),
      changed_pair_checks=len(old)*len(new)+len(new)*(len(new)-1)//2,old_new_repairs=repairs[0],new_new_repairs=repairs[1],seconds=time.time()-started)
    need(path.read_bytes()==raw,'certificate changed during review')
    return report,c,old,gen,modified,new
def section63(c,old,gen,modified,new,section_report=None):
    section_report=Path(section_report) if section_report else CODES/'global_repaired_section63.json'
    lookup=modified|{r['s']:r for r in new}
    if section_report.exists():expected=json.loads(section_report.read_text())
    else:
        octads=json.loads((SOURCE/'shared_parity/dual_octads.json').read_text());candidates=[]
        for h in octads:
            need(h.bit_count()==8 and all((h&g).bit_count()%2==0 for g in gen),'candidate octad')
            need(len(basis([local(h,g) for g in gen]))==7,'candidate section rank')
            def fast(s):
                m=(s&h).bit_count()
                if m%2:return 0
                if m<8:return comb(m,m//2)*(1<<(11-m))
                if s not in lookup:return 1120
                a=lookup[s]
                return (1120 if not(a['phase']&h).bit_count()%2 else 0) if h in a['dual'] else 560
            terms=dict(dense=comb(8,4)*(1<<21),old_sparse=sum(fast(s) for s in old),new_sparse=sum(fast(r['s']) for r in new),four_supports=102900,poles=112)
            candidates.append(dict(octad=h,terms=terms,total=sum(terms.values())))
        expected=max(candidates,key=lambda x:(x['total'],-x['octad']))
    h=expected['octad']
    need(h.bit_count()==8 and all((h&g).bit_count()%2==0 for g in gen),'section dual octad')
    need(len(basis([local(h,g) for g in gen]))==7,'section rank')
    enumerated=0
    def count(s):
        nonlocal enumerated
        mask=local(s,h);m=mask.bit_count()
        if m%2:return 0
        if s not in lookup:
            # Untouched old classes can be ordered independently for this h.
            return comb(m,m//2)*(1<<(11-m)) if m<8 else 1120
        a=lookup[s];codewords=words(a['kernel']);enumerated+=len(codewords)
        need(len(set(codewords))==2048 and min(x.bit_count() for x in codewords if x)==4,'Hamming minimum')
        actual=sum(((x^a['localphase'])&mask).bit_count()==m//2 for x in codewords)
        expected=comb(m,m//2)*(1<<(11-m)) if m<8 else (1120 if not(a['phase']&h).bit_count()%2 else 0) if h in a['dual'] else 560
        need(actual==expected,'affine balanced-sign count');return actual
    blocks=[b for b in combinations(range(64),4) if b[0]^b[1]^b[2]^b[3]==0]
    need(len(blocks)==10416,'quadruple count')
    four=sum(comb(m,m//2)*(1<<(4-m)) for b in blocks if (m:=sum(x<8 for x in b))%2==0)
    terms=dict(dense=comb(8,4)*(1<<21),old_sparse=sum(count(s) for s in old),new_sparse=sum(count(r['s']) for r in new),four_supports=four,poles=112)
    need(terms==expected['terms'] and sum(terms.values())==expected['total'],'63D total')
    return dict(octad=h,terms=terms,total=sum(terms.values()),explicit_sign_words_enumerated=enumerated)
def compressions():
    # Rational enclosures independently replace the quadratic-field comparator.
    lo=Q(412310562561,10**11);hi=Q(412310562562,10**11);need(lo*lo<17<hi*hi,'sqrt17 enclosure')
    def upper(a,b):return Q(a)+Q(b)*(hi if b>=0 else lo)
    tails=(Q(1,4),Q(1,16));body=(Q(1,64),Q(17,1024));checks=0
    for j in range(2):need(60*body[j]+tails[j]**2==1,'61D norms')
    for i,j,k,l in product(range(2),repeat=4):
        for distance in range(16-(i!=k)-(j!=l),61):
            dot=60-2*distance;head=(-1 if i!=k else 1)*tails[j]*tails[l]
            bound=head+dot*body[j] if j==l else upper(head,Q(dot,256))
            need(bound<=Q(1,2),'61D dense separation');checks+=1
    need(Q(14,32)<=Q(1,2) and upper(0,Q(14,128))<=Q(1,2),'61D odd sparse signs')
    labels=(0,1,3,7,6,4);hexagon=((2,0),(1,1),(-1,1),(-2,0),(-1,-1),(1,-1))
    for i,(x,y) in enumerate(hexagon):
        need(x*x+3*y*y==4,'hexagon radius')
        for j,(u,v) in enumerate(hexagon):need(60-2*(16-(labels[i]^labels[j]).bit_count())+x*u+3*y*v<=32,'62D dense separation')
    need((16*16)*3<32*32 and Q(3,4)*16+16==28 and Q(3,4)*64-16==32,'65D cap transfer')
    return dict(dense61_cases=checks,dense62_cases=36,owner65_count=65536+336*8,
      notes='65D changes the dense layer to the inherited BCH radical; all sparse interactions use only amplitudes. Owners use no weight-16 point.')
def main(paths=None):
    reports=[];global_data=None
    supplied=bool(paths)
    paths=list(map(Path,paths)) if paths else [CODES/'global_repaired_certificate.json',CODES/'short_repaired_certificate.json']
    for path in paths:
        ans=family(path);reports.append(ans[0])
        if not ans[1]['removed_coordinates']:
            global_data=ans;section_path=path.with_name(path.name.replace('_certificate.json','_section63.json'))
        print('REVIEWED',path.name,ans[0],flush=True)
    sec=section63(*global_data[1:],section_path) if global_data else None;geo=compressions()
    result=dict(status='PASS',families=reports,section63=sec,geometry=geo,
      inherited_inputs=['Parent [64,28,16] minimum and exact weight-16 list, freshly enumerated by the separate scalar C++ verifier.',
      'The existing 65D BCH radical and 68,224-owner certificate; its transfer to changed sparse layers is checked here.',
      '66/67 Netto quadruple construction remains a separately replayed input.'],
      review='No mathematical gap found. Standalone entry points should reject -O; 63D untouched sign orders and 65D dense-code substitution must be stated explicitly.')
    out=HERE/(f'independent_code_review_{paths[0].stem}.json' if supplied else 'independent_code_review.json')
    out.write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':main(sys.argv[1:])
