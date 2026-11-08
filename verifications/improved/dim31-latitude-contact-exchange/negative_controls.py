"""Reject duplicates, wrong norms and an exact forbidden angle."""
from pathlib import Path
import json,tempfile
import numpy as np
from verify_geometry import verify,require

if __name__=='__main__':
    a=np.zeros(31,np.int16);a[:2]=8
    b=np.zeros(31,np.int16);b[0]=8;b[1:5]=4
    cases={'duplicate':np.array([a,a]),'wrong_norm':np.array([a,a+np.eye(1,31,30,dtype=np.int16)[0]]),'forbidden_angle':np.array([a,b])}
    with tempfile.TemporaryDirectory() as directory:
        for name,points in cases.items():
            path=Path(directory)/(name+'.npz');np.savez_compressed(path,points=points,squared_norm=128,inner_product_limit=64,coordinate_denominator=2)
            for precision in (32,64):
                try:verify(path,precision)
                except ValueError:pass
                else:raise ValueError('Accepted negative control: '+name)
    print(json.dumps({'status':'PASS','rejections':2*len(cases)}))
