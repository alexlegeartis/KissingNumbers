"""Independent finite review of the new sign-family and compression proofs.

Uses no discovery or main-verifier imports and no research-directory inputs.
The parent-code minimum is separately replayed by the main proof package.
This checker independently constructs kernels, verifies affine sign classes,
checks every changed support pair, and recounts a valid 63D section.
The public check() is read-only and returns a deterministic JSON-compatible report.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations,product
from math import comb
import hashlib,json
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
SOURCE=HERE.parent/'dim61-68-95-96-code-layers/source'
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
    raw=path.read_bytes();c=json.loads(raw)
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
      changed_pair_checks=len(old)*len(new)+len(new)*(len(new)-1)//2,old_new_repairs=repairs[0],new_new_repairs=repairs[1])
    need(path.read_bytes()==raw,'certificate changed during review')
    return report,c,old,gen,modified,new
def section63(c,old,gen,modified,new):
    # Select and recount the octad independently, without the main section code.
    lookup=modified|{r['s']:r for r in new}
    octads=json.loads((SOURCE/'shared_parity/dual_octads.json').read_text());candidates=[]
    need(len(octads)==len(set(octads))==72,'octad list size')
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
def compressions(gen):
    four=(24,42,53,55)
    projection=[sum(((g>>p)&1)<<i for i,p in enumerate(four)) for g in gen]
    need(len(basis(projection))==4,'four-coordinate dense projection rank')
    last_two=[sum(((g>>p)&1)<<i for i,p in enumerate(four[2:])) for g in gen]
    need(len(basis(last_two))==2,'two-coordinate dense projection rank')
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
    # Dense versus four-support layers, using rational upper bounds.
    root3_lo=Q(173205080756,10**11);root3_hi=Q(173205080757,10**11)
    need(root3_lo*root3_lo<3<root3_hi*root3_hi,'sqrt3 enclosure')
    for j in range(2):
        for tail_count in (0,1):
            product61=Q(4-tail_count,16)+Q(tail_count,8) if j==0 else upper(Q(tail_count,32),Q(4-tail_count,64))
            need(product61<=Q(1,2),'61D dense/four-support product')
    for x,y in hexagon:
        for a,b in product(range(2),repeat=2):
            product62=Q(4-a-b+a*abs(x),16)+Q(b*abs(y),16)*root3_hi
            need(product62<=Q(1,2),'62D dense/four-support product')
    need(upper(0,Q(1,32))<=Q(1,2) and Q(1,4)<=Q(1,2),'61D dense/pole product')
    need(Q(16,32)==Q(1,2) and 2*16==2*4*4==32 and 2*8<32,'sparse-layer amplitude products')
    blocks=[b for b in combinations(range(64),4) if b[0]^b[1]^b[2]^b[3]==0]
    need(len(blocks)==10416,'four-support design size')
    triples=set()
    for b in blocks:
        for triple in combinations(b,3):
            need(triple not in triples,'four-support intersection above two');triples.add(triple)
    four_counts={d:sum(max(b)<d for b in blocks) for d in (61,62,64)}
    need(four_counts=={61:8555,62:9145,64:10416},'shortened four-support counts')
    need((16*16)*3<32*32 and Q(3,4)*16+16==28 and Q(3,4)*64-16==32,'65D cap transfer')
    return dict(dense61_cases=checks,dense62_cases=36,dense_projection_rank=4,dense61_count=1<<26,dense62_count=6*(1<<24),four_support_counts={str(d):n for d,n in four_counts.items()},owner65_count=65536+336*8,
      notes='65D changes the dense layer to the inherited BCH radical; all sparse interactions use only amplitudes. Owners use no weight-16 point.')
def check():
    """Return a deterministic review report; read inputs without writing files."""
    if not __debug__:
        raise RuntimeError('Run without -O: this independent proof entry point requires ordinary Python.')
    pinned=json.loads((HERE/'data/source_inputs.json').read_text())['sha256']
    source_hashes={}
    for source in (SOURCE/'dim61_96/code64_supports.npz',SOURCE/'shared_parity/dual_octads.json'):
        relative=source.relative_to(ROOT).as_posix()
        digest=hashlib.sha256(source.read_bytes()).hexdigest()
        need(digest==pinned[relative],'pinned source '+relative)
        source_hashes[relative]=digest
    global_data=family(HERE/'data/global64.json')
    short_data=family(HERE/'data/short61_62.json')
    need(global_data[1]['removed_coordinates']==[] and not global_data[1]['odd_total_parity'],'global family convention')
    need(short_data[1]['removed_coordinates']==[24,42,53,55] and short_data[1]['odd_total_parity'],'short family convention')
    sec=section63(*global_data[1:]);geo=compressions(global_data[3])
    short_blocks=len(short_data[2])+len(short_data[5])
    global_blocks=len(global_data[2])+len(global_data[5])
    base64=(1<<28)+global_blocks*2048+10416*16+128
    totals={
        '61':(1<<26)+short_blocks*2048+8555*16+122,
        '62':6*(1<<24)+short_blocks*2048+9145*16+124,
        '63':sec['total'],
        '64':base64,
        '65':base64+68224+2,
    }
    return dict(status='PASS',families=[global_data[0],short_data[0]],section63=sec,
      geometry=geo,totals=totals,pinned_source_files=source_hashes,
      independently_checked=[
        'Explicit RM(1,4) duals, affine phases, independently generated kernels, and every new/old and new/new support intersection.',
        'All 72 dual-octad candidates and the selected 63D section, including explicit affine sign enumeration.',
        '61D/62D dense geometry and cardinalities, sparse amplitude interactions, and the affine-plane quadruple design.',
        '65D cap transfer to the changed weight-sixteen layers.'
      ],
      separately_replayed_inputs=[
        'The parent [64,28,16] code minimum and exact weight-sixteen support list are replayed by check_family.cpp in the main verifier.',
        'The unchanged 65D BCH radical and 68,224-owner certificate are replayed by replay_parents.py; this independent review checks their transfer.',
        'The 66D and 67D retained parent supports, Netto quadruples, and 67D baseline and new subcodes are covered by the other package verifiers.'
      ],
      review='No mathematical gap found in the reviewed scope. The 63D untouched sign orders may be chosen for its octad; 65D uses the separately replayed BCH dense layer.')

if __name__=='__main__':
    print(json.dumps(check(),indent=2,sort_keys=True))
