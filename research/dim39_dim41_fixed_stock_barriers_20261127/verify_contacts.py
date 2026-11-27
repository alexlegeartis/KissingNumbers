from pathlib import Path
import json
import numpy as np
P=Path(__file__).resolve().parent
D=[np.load(P/f"configuration_{i}.npz",allow_pickle=False) for i in (41,42)]
G=D[0]["gram"].astype("int64")
assert np.array_equal(G,D[1]["gram"])
rows=json.loads((P/"contacts41.json").read_text())
assert len(rows)==3
for rec in rows:
    a,b=map(lambda k:np.array(rec[k],dtype="int64"),["X","Y"])
    i,j=rec["xpin"]; k,l=rec["ypin"]
    assert np.array_equal(a%3,D[i]["residues"][j])
    assert np.array_equal(b%3,D[k]["residues"][l])
    assert a@G@a==36 and b@G@b==42 and a@G@b==21
    print("Certified body-product witness",rec["pair"])
