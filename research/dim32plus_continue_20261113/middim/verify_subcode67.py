"""Exact independent verification of the affine 67D sign-subcode extension.

All changed sign words and every affected support interaction are replayed.
The original 67D construction and Netto quadruples are inherited from their
existing verified package; no original point count is increased by assertion.
"""
from pathlib import Path
import hashlib,json,sys
import numpy as np
HERE=Path(__file__).parent;ROOT=HERE.parents[2]
SOURCE=ROOT/'verifications/improved/dim61-68-95-96-code-layers/source'
def need(ok,message):
    if not ok:raise ValueError(message)
def span(rows):
    out=[0]
    for h in rows:out += [x^h for x in out]
    return out
def expand(s,w):
    out=0;j=0
    for i in range(s.bit_length()):
        if s>>i&1:out|=((w>>j)&1)<<i;j+=1
    return out
def kernel(checks):
    piv={};out=[]
    for j in range(16):
        x=sum(((h>>j)&1)<<i for i,h in enumerate(checks));v=1<<j
        while x:
            k=x.bit_length()-1
            if k in piv:x^=piv[k][0];v^=piv[k][1]
            else:piv[k]=(x,v);break
        if not x:out.append(v)
    return out
def code(row,full=False):
    s=row['support'];checks=row['checks'];shift=row['shift']
    need(type(s)is int and s>0 and s.bit_count()==16,'support weight')
    need(0<=shift<65536 and all(type(h)is int and 0<=h<65536 for h in checks),'code range')
    need(all(h.bit_count()%2==0 for h in checks),'antipodality')
    dual=set(span(checks));ker=kernel(checks);linear=span(ker)
    need(len(dual)*len(linear)==65536 and len(linear)==len(set(linear)),'rank identity')
    need(all((w&h).bit_count()%2==0 for w in linear for h in checks),'kernel equations')
    minimum=min(w.bit_count() for w in linear if w)
    need(minimum>=4 and 65535 in linear,'sign minimum and antipodality')
    if full:
        need(len(linear)==2048 and len(checks)==5 and checks[0]==65535,'full Hamming size')
        cols=[sum(((h>>j)&1)<<i for i,h in enumerate(checks)) for j in range(16)]
        need(set(cols)==set(range(1,32,2)),'full Hamming columns')
    elif 'count' in row:need(len(linear)==row['count']==1<<row['dimension'],'new count')
    return dict(s=s,phase=expand(s,shift),dual={expand(s,h) for h in dual},
      count=len(linear),minimum=minimum,words=[expand(s,w^shift) for w in linear],checks=checks,shift=shift)
def load_parent():
    folder=SOURCE/'dim61_96';raw=np.fromfile(folder/'cyclic73.weight16.bin',dtype='<u8').reshape(-1,2)
    choice=json.loads((folder/'shortening_support_search.json').read_text())['dimensions']['67']
    deleted=choice['deleted_coordinates'];dual9=(21,37,38,49,57,63,67,70,72)
    head=[i for i in range(1,73) if i not in dual9]+[0];tail=[i for i in dual9 if i not in deleted];coords=head+tail
    need(len(coords)==len(set(coords))==67,'coordinate dimension')
    selected=[]
    for lo,hi in raw:
        w=int(lo)|(int(hi)<<64)
        if any(w>>i&1 for i in deleted):continue
        selected.append(sum(((w>>c)&1)<<j for j,c in enumerate(coords)))
    need(len(selected)==len(set(selected))==40428 and all(w.bit_count()==16 for w in selected),'parent support count')
    with np.load(folder/'code64_supports.npz',allow_pickle=False) as z:old=set(map(int,z['words']))
    need({w for w in selected if not w>>64}==old and len(old)==30828,'64-coordinate embedding')
    return selected
def verify(path):
    raw=Path(path).read_bytes();c=json.loads(raw)
    need(c['dimension']==67 and c['format']=='affine-sign-subcode-extension-v1','certificate format')
    for name,sha in c['source_sha256'].items():need(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha,'source hash '+name)
    parent=load_parent();parent_set=set(parent)
    fixed_old={};fixed_new=[]
    with np.load(SOURCE/c['base_sign_family'],allow_pickle=False) as z:
        for kind in ('old','new'):
            for s,h,t in zip(z[kind+'_supports'],z[kind+'_duals'],z[kind+'_shifts']):
                row=dict(support=int(s),checks=list(map(int,h)),shift=int(t))
                # Published classes are inherited; rebuild their dual/phase here.
                data=dict(s=int(s),phase=expand(int(s),int(t)),dual={expand(int(s),v) for v in span(list(map(int,h)))},checks=list(map(int,h)),shift=int(t))
                if kind=='old':fixed_old[int(s)]=data
                else:fixed_new.append(data)
    need(len(fixed_old)==len(fixed_new)==1792,'published family size')
    modified={r['support']:code(r,True) for r in c['modified_parent_classes']}
    need(len(modified)==len(c['modified_parent_classes']) and set(modified)<=parent_set,'modified identity')
    need(not set(modified)&set(fixed_old),'previous repairs remain fixed')
    added=[code(row) for row in c['new_classes']];supports={a['s'] for a in added}
    need(len(supports)==len(added) and all(s<1<<64 for s in supports),'new class identity')
    need(not supports&parent_set and not supports&{a['s'] for a in fixed_new},'new supports are absent')
    changed_old_pairs=0
    for s in modified:
        for t in parent:
            if s!=t:need((s&t).bit_count()<=8,'changed old/parent overlap');changed_old_pairs+=s!=t
        for b in fixed_new:need((s&b['s']).bit_count()<=8,'changed old/published new overlap');changed_old_pairs+=1
    old_data=fixed_old|modified;repaired=[];new_pair_checks=0
    def pair(a,b,label):
        nonlocal new_pair_checks
        w=(a['s']&b['s']).bit_count();new_pair_checks+=1
        if w<=8:return
        need(w<=10 and 'dual'in b,'overlap lacks a specified repair')
        intersection=a['s']&b['s']
        common=[h for h in a['dual']&b['dual'] if h.bit_count()==8 and h&intersection==h and ((a['phase']^b['phase'])&h).bit_count()%2]
        need(bool(common),'opposite common eight-check')
        repaired.append(dict(new_support=a['s'],old_support=b['s'],intersection=w,check=min(common),kind=label))
    for i,a in enumerate(added):
        for s in parent:pair(a,old_data.get(s,dict(s=s)),'parent')
        for b in fixed_new:pair(a,b,'published_new')
        for b in added[:i]:pair(a,b,'added')
    terms=json.loads((SOURCE/'continuation49_72/netto.json').read_text())['dimensions']['67']['terms']
    need(terms['dense64']==1<<28 and terms['weight16']==40428*2048 and terms['poles']==134,'base terms')
    base=sum(terms.values())+1792*2048;gain=sum(a['count'] for a in added)
    need(base==c['base_count']==355090822 and base+gain==c['total'],'total')
    # Every new vector has 16 entries of magnitude2: norm64. Other layer
    # products are <=32 (dense), <=32 (four-support), and <=16 (axes).
    need(16*4==64 and 16*2==32 and 4*2*4==32 and 2*8==16,'amplitude bounds')
    return dict(status='PASS',dimension=67,base=base,gain=gain,total=base+gain,
      new_classes=len(added),class_sizes=[a['count'] for a in added],minimum_sign_distances=[a['minimum'] for a in added],
      modified_full_hamming_classes=len(modified),changed_old_support_pairs=changed_old_pairs,
      new_support_pair_checks=new_pair_checks,ten_or_nine_overlap_repairs=repaired,
      certificate_sha256=hashlib.sha256(raw).hexdigest(),
      scope='All changed sign codes and support interactions are freshly exact. The unchanged published 67D family and Netto quadruple construction remain inherited.')
if __name__=='__main__':print(json.dumps(verify(Path(sys.argv[1]) if len(sys.argv)>1 else HERE/'subcode67_certificate.json'),indent=2))
