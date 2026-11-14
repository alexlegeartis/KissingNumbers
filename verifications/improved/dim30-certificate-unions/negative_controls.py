"""Small semantic rejection controls for the production geometry checker."""
from pathlib import Path
import subprocess,sys,tempfile
import numpy as np
ROOT=Path(__file__).resolve().parent
def main():
    A=np.zeros((2,26),np.int16);B=np.zeros_like(A);A[0,[0,1]]=8;A[1,[0,2]]=8
    cases=[]
    cases.append(('valid',A.copy(),B.copy(),3,64,True))
    C=A.copy();C[1]=C[0];cases.append(('duplicate',C,B.copy(),3,64,False))
    C=A.copy();C[0,0]+=1;cases.append(('wrong norm',C,B.copy(),3,64,False))
    C=A.copy();C[1]=0;C[1,0]=8;C[1,[1,2,3,4]]=4
    cases.append(('distinct rational conflict',C,B.copy(),3,64,False))
    C=A.copy();D=B.copy();C[1]=0;C[1,[24,25]]=4;D[1,[0,1]]=4
    cases.append(('radical conflict',C,D,3,64,False))
    cases.append(('wrong field',A.copy(),B.copy(),2,64,False))
    cases.append(('relaxed threshold',A.copy(),B.copy(),3,65,False))
    with tempfile.TemporaryDirectory() as td:
        for i,(name,C,D,r,limit,ok) in enumerate(cases):
            p=Path(td)/f'{i}.npz';o=Path(td)/f'{i}.json'
            np.savez(p,rational=C,radical=D,radicand=r,squared_norm=128,
                     inner_product_limit=limit,coordinate_denominator=2)
            proc=subprocess.run([sys.executable,str(ROOT/'verify_geometry.py'),str(p),'--out',str(o)],capture_output=True,text=True)
            assert (proc.returncode==0)==ok,(name,proc.stdout,proc.stderr)
            print(('ACCEPTED' if ok else 'REJECTED'),name,flush=True)
    print('PASS: one valid control and six semantic corruptions.')
if __name__=='__main__':main()
