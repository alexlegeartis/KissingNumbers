from pathlib import Path
import json,hashlib
import numpy as np
from module_recipe import reconstruct,canonical_hash
from verify_geometry import verify
D=Path(__file__).resolve().parent/'data'
def main():
 recipe=json.load(open(D/'construction.json'));m,_,_,_=reconstruct(recipe)
 z=np.load(D/'module.npz');T=z['transformation_numerator'];d=int(z['transformation_denominator']);assert np.array_equal(T@T.T,d*d*np.eye(15,dtype=np.int64));assert not np.any((m@T.T)%d)
 y=m@T.T//d;assert np.array_equal(np.unique(y,axis=0),np.unique(z['points'],axis=0))
 source=np.load(D/'source31_251395.npz')['points'].astype(np.int64)//2;final=np.load(D/'coordinates31_251396.npz')['points'].astype(np.int64)//2;ex=np.load(D/'exchange_norm32.npz');added=ex['added'].astype(np.int64);removed=ex['removed'].astype(np.int64)
 keys=lambda a:{tuple(map(int,row)) for row in a}
 sk,fk,ak,rk=map(keys,[source,final,added,removed]);assert len(sk)==251395 and len(fk)==251396 and len(ak)==1004 and len(rk)==1003
 assert rk<=sk and not (ak&sk) and fk==(sk-rk)|ak
 E=[0,2,4,5,6,10,11,23]+list(range(24,31));outside=[i for i in range(31) if i not in E];assert not np.any(added[:,outside]);assert keys(added[:,E])<=keys(y)
 assert np.all(np.sum(final**2,axis=1)==32)
 finalindex={tuple(row):i for i,row in enumerate(final)}
 for row in added:
  dots=final@row;dots[finalindex[tuple(row)]]=-32;assert dots.max()<=16
 report=verify(D/'coordinates31_251396.npz',precision=64);report.update(exact_changed_pair_method='direct int64',added=1004,removed=1003,complete_constructed_module=2564)
 (D/'verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':main()
