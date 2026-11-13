"""Four tiny malformed witnesses must be rejected by the geometry verifier."""
import tempfile,json
from pathlib import Path
import numpy as np
from verify_geometry import verify

def main():
    a=np.zeros((2,29),dtype=np.int16);a[0,:2]=8;a[1,2:4]=8
    variants={}
    b=a.copy();b[0,0]=7;variants['wrong_norm']=(b,128)
    b=a.copy();b[1]=b[0];variants['duplicate']=(b,128)
    b=a.copy();b[1]=0;b[1,:5]=[10,4,2,2,2];variants['illegal_angle']=(b,128)
    variants['wrong_metadata']=(a,127)
    rejected=[]
    with tempfile.TemporaryDirectory() as td:
        for name,(x,norm) in variants.items():
            f=Path(td)/'bad.npz';np.savez(f,points=x,squared_norm=norm,inner_product_limit=64,coordinate_denominator=2)
            try:verify(f,precision=64)
            except ValueError:rejected.append(name)
            else:raise RuntimeError('Failed to reject '+name)
    print(json.dumps({'pass':True,'rejected':rejected}))
if __name__=='__main__':main()
