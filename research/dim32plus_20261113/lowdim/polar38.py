"""Exact enumeration of sparse 38D latitude heads; no kissing-bound claim."""
from pathlib import Path
from itertools import combinations
import json
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
P=ROOT/'verifications/improved/dim38-joint-latitudes'
W=np.loadtxt(P/'source/source/dim14_directions.txt',delimiter=',',dtype=np.int16)
Q=np.array(json.loads((P/'source/certificate.json').read_text())['rotation_numerator'],np.int16)
F=W[np.count_nonzero(W,axis=1)==8]@Q.T
signs=1-2*((np.arange(32)[:,None]>>np.arange(5))&1)
signs=np.column_stack([np.ones(32,dtype=int),signs]).astype(np.int16)
accepted=[];hist={};tested=0
for support in combinations(range(14),6):
 Y=np.zeros((32,14),np.int16);Y[:,support]=signs
 maxcap=np.abs(Y@W.T).max(axis=1)
 maxaxis=np.abs(Y@F.T).max(axis=1)
 valid=(maxcap<=4)&(maxaxis<=4)
 for i in np.flatnonzero(valid):accepted.append(Y[i])
 if valid.any():hist[','.join(map(str,support))]=int(valid.sum())
 tested+=32
np.savez_compressed(Path(__file__).with_name('polar38_heads.npz'),heads=np.array(accepted))
report={'oriented_sign_lines_tested':tested,'accepted_head_lines':len(accepted),'tail_support_count':len(hist),'sign_counts_per_support':hist,'scope':'Complete weight-six ±1 head screen with sufficient cap and fixed-flat-axis bound |dot|≤4; no body allocation or kissing improvement follows.'}
Path(__file__).with_name('polar38_report.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='sign_counts_per_support'}))
