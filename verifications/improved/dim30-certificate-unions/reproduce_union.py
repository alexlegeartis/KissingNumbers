"""Optional discovery replay; requires SciPy as well as NumPy."""
from pathlib import Path
import numpy as np
from union_engine import read,union
ROOT=Path(__file__).resolve().parent

def main():
    A,B,k=read(ROOT/'data/source30_248180.npz')
    for n,expected in ((248020,248204),(247972,248212)):
        C,D,kk=read(ROOT/f'data/donor30_{n}.npz')
        assert not np.any(B) and not np.any(D)
        A,B,report=union(A,B,C,D,3)
        assert len(A)==expected
        print(n,report,flush=True)
    frozen=np.load(ROOT/'data/coordinates30_248212.npz')['rational']
    assert {v.tobytes() for v in A}=={v.tobytes() for v in frozen}
    print('PASS: both exact matching exchanges reproduce the frozen code.')
if __name__=='__main__':main()
