"""Independent direct numerical pair replay of saved wave point sets."""
from pathlib import Path
import hashlib,json,time
import numpy as np
P=Path(__file__).resolve().parent

def main():
    started=time.monotonic();rows=[]
    for path in sorted(P.glob('*.npz')):
        data=np.load(path)
        if 'vectors' not in data:continue
        X=data['vectors']
        if not np.isfinite(X).all():raise ValueError('Nonfinite '+path.name)
        maximum=-1.;violations=0;duplicates=0
        for i in range(len(X)):
            dots=X[i+1:]@X[i]
            if len(dots):
                maximum=max(maximum,float(dots.max()));violations+=int(np.count_nonzero(dots>.5+1e-10));duplicates+=int(np.count_nonzero(dots>1-1e-10))
        normerr=float(np.max(abs(np.sum(X*X,axis=1)-1)))
        rows.append({'file':path.name,'points':len(X),'dimension':X.shape[1],'pairs_checked':len(X)*(len(X)-1)//2,'maximum_cosine':maximum,'violating_pairs_at_1e-10':violations,'near_duplicate_pairs':duplicates,'norm_squared_error':normerr,'numerically_admissible':maximum<=.5+1e-10 and normerr<1e-10,'exact_certificate':False,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    result={'method':'separate direct row-vs-later-rows dot products; no orbit compression or objective reuse','point_sets':rows,'seconds':time.monotonic()-started,'status':'Numerical screening only; no exact or rigorous interval certificate.'}
    (P/'replay_results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
