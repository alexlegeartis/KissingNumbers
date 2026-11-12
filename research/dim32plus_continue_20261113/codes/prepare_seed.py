"""Convert the published sign certificates to exact search seeds."""
from pathlib import Path
import itertools
import json
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parents[2]/'verifications/improved/dim61-68-95-96-code-layers/source'
sys.path.insert(0,str(SOURCE/'shared_parity'))
from shared_layers import words,expand,solve

DTYPE=np.dtype([('a','<u8'),('old','<u4',(8,))])

def shape(a,indices,old):
    E=[a&~int(old[j]) for j in indices]
    L=[j for j in range(8) if E[0]&E[j]]
    R=[j for j in range(8) if not E[0]&E[j]]
    assert len(L)==len(R)==4
    O=0
    for j in L:O|=E[j]
    Q=a^O
    L=sorted(L,key=lambda j:O^E[j]);R=sorted(R,key=lambda j:Q^E[j])
    P=[O^E[j] for j in L];T=[Q^E[j] for j in R]
    def checks(pm):
        H=[0]*8
        for j in range(4):
            H[L[j]]=P[j]|(Q^T[pm[j]])
            H[R[pm[j]]]=T[pm[j]]|(O^P[j])
        return H
    return checks

def main():
    old=np.fromfile(HERE/'old_supports.bin',dtype='<u8')
    records=np.fromfile(HERE/'candidates.bin',dtype=DTYPE)
    with np.load(SOURCE/'shared_parity/global.npz',allow_pickle=False) as z:
        data={}
        for kind in ('old','new'):
            for s,d,t in zip(z[kind+'_supports'],z[kind+'_duals'],z[kind+'_shifts']):
                s=int(s);d=list(map(int,d))
                data[s]=dict(checks=d,span={expand(s,h):i for i,h in enumerate(words(d))},phase=expand(s,int(t)))
        seeds=list(map(int,z['new_supports']))
    masks={'global64':0,'short61':sum(1<<c for c in (24,42,53,55))}
    report={}
    for name,removed in masks.items():
        if len(sys.argv)>1 and name!=sys.argv[1]:continue
        pool=records[(records['a']&np.uint64(removed))==0]
        pool.tofile(HERE/f'{name}_candidates.bin')
        ids={int(r['a']):i for i,r in enumerate(pool)}
        selected=[s for s in seeds if not s&removed]
        phases={s:x['phase'] for s,x in data.items()}
        if removed:
            supports=sorted(s for s in data if not s&removed)
            pos={s:i for i,s in enumerate(supports)}
            equations=[]
            for a in selected:
                for j in pool[ids[a]]['old']:
                    b=int(old[j])
                    if b&removed:continue
                    common=[h for h in data[a]['span'] if h.bit_count()==8 and h&(a&b)==h and h in data[b]['span']]
                    assert len(common)==1
                    h=common[0];ca=data[a]['span'][h];cb=data[b]['span'][h]
                    equations.append((((ca>>1)<<(4*pos[a]))^((cb>>1)<<(4*pos[b])),1^(ca&1)^(cb&1)))
            value,rank=solve(equations);assert value is not None
            for s in supports:
                target=1|(((value>>(4*pos[s]))&15)<<1)
                phase,_=solve([(h,(target>>j)&1) for j,h in enumerate(data[s]['checks'])])
                assert phase is not None and phase.bit_count()%2==1
                phases[s]=expand(s,phase)
        rows=[]
        for s in selected:
            i=ids[s];checks=shape(s,pool[i]['old'],old)
            pm=next(p for p in itertools.permutations(range(4)) if all(h in data[s]['span'] for h in checks(p)))
            rows.append(' '.join(map(str,(s,phases[s],i,*pm))))
        (HERE/f'{name}_best.txt').write_text(str(len(rows))+'\n'+'\n'.join(rows)+'\n')
        report[name]=dict(candidates=len(pool),seed_supports=len(rows),removed_mask=removed,odd=bool(removed))
    if len(sys.argv)==1:
        (HERE/'seed_preparation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)

if __name__=='__main__':main()
