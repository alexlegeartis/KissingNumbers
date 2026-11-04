"""Weight exact conflicts by their actual sign-orbit squared overlap energy."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import hashlib,json
import numpy as np
from numba import njit,prange,set_num_threads
from search import ROOT


@njit(cache=True,parallel=True)
def weighted(bad,reps):
    n=len(reps);result=np.zeros((n,n),dtype=np.uint8)
    for i in prange(n):
        for j in range(i):
            if not bad[i,j]:continue
            dots=np.zeros(3)
            for b in range(3):
                for d in range(4):dots[b]+=reps[i,4*b+d]*reps[j,4*b+d]
            energy=0.
            for s in range(8):
                dot=0.
                for b in range(3):dot+=dots[b]*(1 if s&(1<<b) else -1)
                v=max(dot-.5,0.);energy+=8*v*v
            # Positive integer for every exact conflict; zero iff compatible.
            w=min(255,max(1,int(round(energy*128))))
            result[i,j]=w;result[j,i]=w
    return result


def main():
    set_num_threads(96);out=ROOT/'runs'/'orbit_pool';n=20736
    bad=np.fromfile(out/'conflicts.bin',dtype=np.uint8).reshape(n,n)
    reps=np.fromfile(out/'representatives.bin',dtype=np.float64).reshape(n,12)
    w=weighted(bad,reps)
    assert np.array_equal(w>0,bad.astype(bool))
    w.tofile(out/'weights.bin')
    r=dict(sha256=hashlib.sha256(w.tobytes()).hexdigest(),
           weight='max(1,round(128 * full sign-orbit squared violation)), capped at 255',
           zero_pattern_exactly_matches_conflict_graph=True)
    (out/'weights_report.json').write_text(json.dumps(r,indent=2)+'\n');print(r)


if __name__=='__main__':main()
