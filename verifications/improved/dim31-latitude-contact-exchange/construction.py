"""Integer block formula for an orthogonal square-to-cube exchange in R31."""
from pathlib import Path
from itertools import product
import hashlib,json
import numpy as np
ROOT=Path(__file__).resolve().parent

def require(ok,why):
    if not ok: raise ValueError(why)

def delta():
    A=np.zeros(31,np.int16);A[[0,2,6,23]]=[-4,4,4,-4]
    B=np.zeros(31,np.int16);B[[8,9,13,14,15,17,18,22]]=[-2,2,2,2,-2,2,2,2]
    C=np.zeros(31,np.int16);C[[24,25]]=[-4,4]
    D=np.zeros(31,np.int16);D[[27,29]]=[4,-4]
    blocks=np.array([A,B,C,D],np.int64)
    require(np.array_equal(blocks@blocks.T,np.diag([64,32,32,32])), 'Orthogonal blocks and norms')
    removed=np.array([a*A+b*(C+D) for a,b in product((1,-1),repeat=2)],np.int16)
    added=np.array([a*A+b*B+c*C for a,b,c in product((1,-1),repeat=3)],np.int16)
    partners=np.array([2*B,-2*B,2*C,-2*C],np.int16)
    return removed,added,partners

def reconstruct():
    row=json.loads((ROOT/'data/provenance.json').read_text())['records'][0]
    source=ROOT/row['baseline'];require(hashlib.sha256(source.read_bytes()).hexdigest()==row['baseline_sha256'],'Pinned baseline bytes')
    old=np.load(source,allow_pickle=False)['points'];removed,added,partners=delta()
    require(old.shape==(251196,31) and old.dtype==np.int16,'Baseline shape and type')
    oldset={v.tobytes() for v in old};rs={v.tobytes() for v in removed};ads={v.tobytes() for v in added}
    require(len(oldset)==len(old) and rs<=oldset and not ads&oldset,'Exact deletion/insertion membership')
    require(all(v.tobytes() in oldset-rs for v in partners),'Four retained A3 partner roots')
    keep=np.array([v.tobytes() not in rs for v in old]);final=np.r_[old[keep],added]
    require(final.shape==(251200,31),'Reconstructed dimension/count')
    return final
