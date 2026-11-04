"""Exact normalized-ray check of every archived 48-system, including axes."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import io,json,zipfile
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent


def check(x):
    vectors=[];norms=[]
    for row in x:
        ratios=[float(t).as_integer_ratio() for t in row]
        exponent=max(d.bit_length()-1 for _,d in ratios)
        z=[n<<(exponent-(d.bit_length()-1)) for n,d in ratios]
        vectors.append(z);norms.append(sum(a*a for a in z))
    coordinate=sum(4*z*z>norms[i] for i,row in enumerate(vectors) for z in row)
    pair=0
    for i,row in enumerate(vectors):
        for j in range(i):
            dot=sum(a*b for a,b in zip(row,vectors[j]))
            if dot>0 and 4*dot*dot>norms[i]*norms[j]:pair+=1
    return dict(passed=coordinate==0 and pair==0,coordinate_violations=int(coordinate),pair_violations=pair)


def main():
    out=ROOT/'runs'/'archive_valid';out.mkdir(parents=True,exist_ok=True);rows=[];good=[]
    with zipfile.ZipFile(ROOT/'inputs'/'saved_48_systems.zip') as z:
        for name in z.namelist():
            if not name.endswith('.npy'):continue
            raw=np.load(io.BytesIO(z.read(name)),allow_pickle=False).astype(np.float64)
            result=dict(name=name,**check(raw));rows.append(result)
            if result['passed']:
                np.save(out/f'valid_{len(good):03d}.npy',raw);good.append(result)
    report=dict(arrays=len(rows),exactly_valid_with_signed_axes=len(good),valid=good,results=rows)
    (ROOT/'archive_exact.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='results'},indent=2))


if __name__=='__main__':main()
