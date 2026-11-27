from pathlib import Path
import json
import numpy as np
p=Path(__file__).resolve().parent
data=[]
fields="gram residues radii heads denominator populations".split()
for dimension in (41,42):
    with np.load(p/f"configuration_{dimension}.npz",allow_pickle=False) as D:
        data.append({k: {"a":D[k].astype(np.int64).ravel().tolist()} for k in fields})
print(json.dumps(data))
