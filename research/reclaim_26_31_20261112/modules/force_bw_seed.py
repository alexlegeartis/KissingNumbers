"""Materialize and exactly check a forced BW15 whole-module replacement."""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
from bw15_transplant import NORMAL

def main():
    ap=argparse.ArgumentParser();ap.add_argument('source',type=Path);ap.add_argument('module',type=Path)
    ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    # Copy module first so a concurrently improving discovery cannot change it.
    raw=a.module.read_bytes();(a.out/'module.npz').write_bytes(raw)
    z=np.load(a.source);X=z['points'].astype(np.int16)
    if int(z['squared_norm'])==128:X//=2
    assert X.shape[1]==31 and np.all(np.sum(X.astype(np.int64)**2,axis=1)==32)
    head=np.load(a.out/'module.npz')['heads'].astype(np.int16)
    assert head.shape==(2016,15)
    Y=np.zeros((len(head),31),np.int16);Y[:,NORMAL]=head[:,:8];Y[:,24:]=head[:,8:]
    assert np.all(np.sum(Y.astype(np.int64)**2,axis=1)==32)
    q=np.sum(X[:,24:].astype(np.int64)**2,axis=1)
    module_old=(q==16)&~np.any(X[:,[j for j in range(24) if j not in NORMAL]],axis=1)
    oldbad=np.zeros(len(X),bool);xf=X.astype(np.float32)
    assert 31*int(np.abs(Y).max())*int(np.abs(X).max())<2**24
    for at in range(0,len(Y),64):oldbad|=np.any(Y[at:at+64].astype(np.float32)@xf.T>16,axis=0)
    D=Y.astype(np.int32)@Y.astype(np.int32).T;np.fill_diagonal(D,-32)
    assert D.max()<=16
    removed_mask=module_old|oldbad;F=np.r_[X[~removed_mask],Y]
    assert len(np.unique(F,axis=0))==len(F)
    source_set={r.tobytes() for r in X};final_set={r.tobytes() for r in F}
    removed=np.asarray([r for r in X if r.tobytes() not in final_set],np.int16)
    added=np.asarray([r for r in F if r.tobytes() not in source_set],np.int16)
    assert len(F)==len(X)+len(added)-len(removed)
    np.savez_compressed(a.out/f'forced31_{len(F)}.npz',points=2*F,squared_norm=128,inner_product_limit=64,coordinate_denominator=2)
    np.savez_compressed(a.out/'exchange.npz',removed=2*removed,added=2*added,coordinate_denominator=2,squared_norm=128,inner_product_limit=64)
    report={'source':str(a.source.resolve()),'source_sha256':hashlib.sha256(a.source.read_bytes()).hexdigest(),
        'module_sha256':hashlib.sha256(raw).hexdigest(),'baseline':len(X),'forced_seed':len(F),'gain':len(F)-len(X),
        'removed':len(removed),'added':len(added),'old_module_heads':int(module_old.sum()),
        'external_deleted':int(np.sum(oldbad&~module_old)),'checks':{'all_source_new_pairs':len(X)*len(Y),'all_new_new_pairs':len(Y)*(len(Y)-1)//2,
            'arithmetic':'bounded exact integer binary32 source-new; independent int32 new-new','new_new_maximum':int(D.max())},
        'scope':'All changed inequalities checked exactly; unchanged source validity must be independently established. This may be below the incumbent and is a research seed.'}
    (a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':main()
