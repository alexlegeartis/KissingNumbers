"""Independent exact replay of the full-row support obstruction, not optimality."""
from pathlib import Path
import hashlib,itertools,json,math
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'verifications/improved/dim88-96-z8-golay/data/supports.npz'
HERE=Path(__file__).resolve().parent
EXPECTED_SOURCE='d0def7d6587ca3aa1b39c9f3e7ccb844310285a7f1794dc8e9492f0fff9fc65d'
EXPECTED_PROOF='6c510520d410b67fa82205f7b018727a2ea56c5fb3cb490e687c091f8b7a8260'

def check(proof=HERE/'full_row_witnesses.bin'):
    if hashlib.sha256(BASE.read_bytes()).hexdigest()!=EXPECTED_SOURCE: raise ValueError('Source hash mismatch')
    if hashlib.sha256(Path(proof).read_bytes()).hexdigest()!=EXPECTED_PROOF: raise ValueError('Proof hash mismatch')
    z=np.load(BASE)
    lo,hi=z['supports_low'],z['supports_high']
    if lo.shape!=(26022,) or hi.shape!=(26022,): raise ValueError('Unexpected support family')
    if lo.dtype!=np.uint64 or hi.dtype!=np.uint64 or np.any(hi>=2**32): raise ValueError('Invalid support encoding')
    if len(np.unique(np.column_stack([lo,hi]),axis=0))!=26022: raise ValueError('Duplicate source supports')
    rows=np.column_stack([((lo>>np.uint64(c))&1) for c in range(64)]+[((hi>>np.uint64(c))&1) for c in range(32)]).reshape(-1,24,4).sum(2)
    if np.any(rows.sum(1)!=24): raise ValueError('Source support weight is not24')
    w=np.fromfile(proof,dtype='<u4')
    if w.shape!=(math.comb(24,6),) or np.any(w>=len(lo)): raise ValueError('Incomplete or invalid witness table')
    count=0; smallest=24
    for index,subset in enumerate(itertools.combinations(range(24),6)):
        overlap=sum(int(rows[w[index],j]) for j in subset)
        if overlap<=12: raise ValueError(f'Invalid blocker at {index}')
        smallest=min(smallest,overlap); count+=1
    return {'status':'exact_bounded_obstruction','candidate_family':'all six full rows of four cells','tested':count,'minimum_blocker_intersection':smallest,'existing_support_sha256':hashlib.sha256(BASE.read_bytes()).hexdigest(),'proof_sha256':hashlib.sha256(Path(proof).read_bytes()).hexdigest(),'scope':'Every candidate violates the sufficient support-intersection<=12 criterion with an existing support. This does not itself prove a signed spherical incompatibility; sign restrictions, replacements and other constructions remain untested.'}
if __name__=='__main__': print(json.dumps(check(),indent=2))
