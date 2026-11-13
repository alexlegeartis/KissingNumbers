"""Bounded joint insertion with complete old Hamming-code preservation."""
from pathlib import Path
import json,random,time,math
from utils import expand,span,need
HERE=Path(__file__).parent

def local(s,w):
    out=0;j=0
    for i in range(s.bit_length()):
        if s>>i&1:out|=((w>>i)&1)<<j;j+=1
    return out
def solve(equations):
    piv={}
    for x,y in equations:
        while x:
            j=x.bit_length()-1
            if j in piv:x^=piv[j][0];y^=piv[j][1]
            else:piv[j]=(x,y);break
        if not x and y:return None
    out=0
    for j in sorted(piv):
        x,y=piv[j]
        if ((x&out).bit_count()%2)!=y:out|=1<<j
    return out
def full_hamming(s,equations):
    shift=solve(equations)
    if shift is None:return None
    rows=[65535];space={0,65535}
    for h,_ in equations:
        if h in space:continue
        rows.append(h);space|={x^h for x in list(space)}
        if len(rows)>5 or any(x.bit_count()not in(0,8,16) for x in space):return None
    groups={}
    for j in range(16):
        key=sum(((h>>j)&1)<<i for i,h in enumerate(rows[1:]));groups.setdefault(key,[]).append(j)
    remaining=5-len(rows)
    if len(groups)!=1<<(len(rows)-1) or any(len(g)!=1<<remaining for g in groups.values()):return None
    for bit in range(remaining):rows.append(sum(1<<j for g in groups.values() for i,j in enumerate(g) if i>>bit&1))
    cols=[sum(((h>>j)&1)<<i for i,h in enumerate(rows)) for j in range(16)]
    need(set(cols)==set(range(1,32,2)),'completion columns')
    need(all((shift&h).bit_count()%2==y for h,y in equations),'completion phase')
    return dict(support=s,checks=rows,shift=shift)
def main(seconds=60):
    frozen=dict(format='affine-sign-subcode-extension-v1',dimension=64,base_count=335711104);context=json.loads((HERE/'context.json').read_text())
    pool=[json.loads(line) for line in (HERE/'feasible.jsonl').read_text().splitlines()]
    need(len({r['support'] for r in pool})==len(pool)==35,'feasible pool')
    for r in pool:
        phase=expand(r['support'],r['phase']);r['requirements']={}
        for kind,index,h in r['free']:
            s=context['old_order'][index] if kind==0 else context['extra_full'][index]
            r['requirements'].setdefault(s,[]).append((local(s,h),1^((phase&h).bit_count()%2)))
    rng=random.Random(20261113);best=0;iteration=0;start=time.time()
    while time.time()-start<seconds:
        noise=[0,.25,.7,1.5,3][iteration%5]
        order=sorted(range(len(pool)),key=lambda i:pool[i]['count']*math.exp(noise*rng.uniform(-1,1)),reverse=True)
        selected=[];constraints={};completions={};value=0
        for i in order:
            r=pool[i]
            if any((r['support']&pool[j]['support']).bit_count()>8 for j in selected):continue
            updates={};good=True
            for s,eq in r['requirements'].items():
                all_eq=constraints.get(s,[])+eq;completion=full_hamming(s,all_eq)
                if completion is None:good=False;break
                updates[s]=(all_eq,completion)
            if not good:continue
            selected.append(i);value+=r['count']
            for s,(eq,comp) in updates.items():constraints[s]=eq;completions[s]=comp
        if value>best:
            best=value;out=dict(frozen);out['new_classes']=[dict(support=pool[i]['support'],checks=pool[i]['checks'],shift=pool[i]['phase'],dimension=pool[i]['dimension'],count=pool[i]['count']) for i in selected]
            out['modified_parent_classes']=[completions[s] for s in sorted(completions)];out['total']=out['base_count']+best
            (HERE/'certificate.json').write_text(json.dumps(out,indent=2)+'\n')
            (HERE/'search.json').write_text(json.dumps(dict(gain=best,classes=len(selected),selected_indices=selected,modified_classes=len(completions),iteration=iteration,seconds=time.time()-start),indent=2))
            print('BEST',best,'classes',len(selected),'old modifications',len(completions),'iteration',iteration,flush=True)
        iteration+=1
    print('DONE',best,'iterations',iteration,'seconds',time.time()-start,flush=True)
if __name__=='__main__':main(15)
