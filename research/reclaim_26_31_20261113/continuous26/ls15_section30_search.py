"""Bounded E8/tail orientation search for 1,356-point LS15 sections in30D."""
from pathlib import Path
import argparse,hashlib,json,math,sys,time
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'packingstar_audit'))
from bw14_search import mutate


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--seconds',type=float,default=120);ap.add_argument('--seed',type=int,default=20261113)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True);start=time.monotonic();root=HERE.parents[2]
    src=root.parent/'PackingStar/25D-31D/new_bounds_20261111/D30_247992_coordinates.npz'
    z=np.load(src);X=z['points'].astype(np.int16);assert int(z['coordinate_denominator'])==1
    M=np.load(HERE.parent/'algebraic/ls15_exact/module.npz')['points'].astype(np.int16)
    J=np.r_[[0,2,4,5,6,10,11,23],np.arange(24,30)];other=[j for j in range(24) if j not in J]
    oldmodule=~np.any(X[:,other],axis=1);proj=X[:,J]
    q=np.sum(X[:,24:].astype(np.int32)**2,axis=1);h=np.sum(X[:,J[:8]].astype(np.int32)**2,axis=1)
    external=(((q==8)&(h>8))|((q==16)&(h>0)))&~oldmodule
    A0=proj[external];assert len(A0)==152
    full0=M[M[:,14]==0][:,:14];assert len(full0)==1356
    headmask=np.sum(full0[:,8:].astype(np.int32)**2,axis=1)==16
    base=full0[headmask];assert len(base)==1056
    assert np.all(np.sum(base[:,:8].astype(np.int32)**2,axis=1)==16)
    basef=base.astype(np.float32);A=A0.copy();R=4*np.eye(14,dtype=np.int16)
    def cost(A):return int(np.any(basef@A.astype(np.float32).T>16,axis=0).sum())
    score=cost(A);best=score;bestA=A.copy();bestR=R.copy();rng=np.random.default_rng(a.seed);records=[];bestgain=-100000
    allprofiles,inv,weights=np.unique(proj,axis=0,return_inverse=True,return_counts=True)
    def record(it):
        nonlocal bestgain
        num=full0.astype(np.int32)@R.astype(np.int32).T;assert np.all(num%4==0);Y=(num//4).astype(np.int16)
        assert np.all(np.sum(Y.astype(np.int32)**2,axis=1)==32)
        blocked=np.any(Y.astype(np.float32)@allprofiles.astype(np.float32).T>16,axis=0);remove=blocked[inv]
        gain=len(Y)-int(remove.sum());ext=int((remove&~oldmodule).sum())
        assert ext==score, 'Compressed exterior price must equal complete source price'
        row=dict(iteration=it,external=ext,removed=int(remove.sum()),added=len(Y),gain=gain,total=len(X)+gain,seconds=time.monotonic()-start)
        if gain>bestgain:
            bestgain=gain;add=np.zeros((len(Y),30),np.int16);add[:,J]=Y;F=np.vstack([X[~remove],add]);assert len(F)==len(np.unique(F,axis=0))
            np.savez_compressed(a.out/'best_module.npz',points=Y,squared_norm=32,inner_product_limit=16,transformation_numerator=R,transformation_denominator=4)
            np.savez_compressed(a.out/f'coordinates30_{len(F)}.npz',points=2*F,squared_norm=128,inner_product_limit=64,coordinate_denominator=2)
            np.savez_compressed(a.out/'exchange.npz',removed=2*X[remove],added=2*add)
        records.append(row);print(json.dumps(row),flush=True)
        (a.out/'report.json').write_text(json.dumps(dict(source=str(src),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),baseline=len(X),old_module=int(oldmodule.sum()),new_module=1356,new_heads=1056,omitted_tail=6,best_gain=bestgain,records=records,scope='Exact source-compressed full-module placement; direct independent changed-pair and full-configuration replay required before record registration'),indent=2)+'\n')
    record(0);it=0
    while time.monotonic()-start<a.seconds:
        it+=1;u=rng.random()
        if u<.5:
            kind='root';choice=rng.choice([-1,1],8);choice[-1]=np.prod(choice[:-1])
        elif u<.65:kind='core_swap';choice=rng.choice(8,2,replace=False)
        elif u<.85:kind='tail_swap';choice=rng.choice(6,2,replace=False)
        else:kind='tail_sign';choice=int(rng.integers(6))
        C=mutate(A,kind,choice);value=cost(C);temp=5*(1-(it%2500)/2500)+.1
        if value<=score or rng.random()<math.exp((score-value)/temp):
            A,R,score=C,mutate(R,kind,choice),value
            if score<best:
                best,bestA,bestR=score,A.copy(),R.copy();record(it)
        if it%2500==0:A,R,score=bestA.copy(),bestR.copy(),best
        if best==0:break
    print(json.dumps(dict(finished=True,iterations=it,best_external=best,best_gain=bestgain,seconds=time.monotonic()-start)),flush=True)


if __name__=='__main__':main()
