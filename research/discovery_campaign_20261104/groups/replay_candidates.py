"""Independent full pair numerical replay of every saved trial point set."""
from pathlib import Path
import json,time,hashlib
import numpy as np
P=Path(__file__).resolve().parent

def main():
    started=time.perf_counter();reports=[]
    for path in sorted(P.glob('*.npz')):
        data=np.load(path)
        if 'vectors' not in data:continue
        X=data['vectors'];maxdot=-1.0
        for i in range(0,len(X),128):
            A=X[i:i+128]@X.T
            A[np.arange(len(A)),i+np.arange(len(A))]=-np.inf
            maxdot=max(maxdot,float(A.max()))
        assert np.isfinite(X).all()
        normerr=float(np.max(abs(np.sum(X*X,axis=1)-1)))
        distinct=len(np.unique(np.round(X,10),axis=0));assert distinct==len(X)
        assert normerr<1e-12
        reports.append({'file':path.name,'dimension':X.shape[1],'points':len(X),'distinct_rounded_10_decimals':distinct,'max_norm_squared_error':normerr,'maximum_pair_cosine':maxdot,'admissible_at_1e-9_numerical_tolerance':maxdot<=.5+1e-9,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'exact_certificate':False})
    result={'replay_method':'Direct full pair dot products from saved vectors, independent of seed objective/orbit compression. Floating screening is not certification.','point_sets':reports,'seconds':time.perf_counter()-started}
    (P/'replay_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
