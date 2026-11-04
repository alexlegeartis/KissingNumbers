"""Certify the actual printed decimal coordinates, without floating arithmetic."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from decimal import Decimal
import hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'inputs'/'gram841_coordinates.txt'


def main():
    start=time.monotonic();vectors=[];norms=[]
    for line in SOURCE.read_text().splitlines():
        tokens=line.split()
        if not tokens:continue
        assert len(tokens)==12
        decimals=[Decimal(t).as_tuple() for t in tokens]
        scale=max(0,max(-t.exponent for t in decimals))
        z=[]
        for t in decimals:
            n=int(''.join(str(i) for i in t.digits))*(-1 if t.sign else 1)
            z.append(n*10**(scale+t.exponent))
        norm=sum(t*t for t in z);assert norm>0;vectors.append(z);norms.append(norm)
    assert len(vectors)==841
    p=2**23-1;q=2**24;violations=0
    for i,row in enumerate(vectors):
        for j in range(i):
            dot=sum(a*b for a,b in zip(row,vectors[j]))
            if dot>0 and q*q*dot*dot>p*p*norms[i]*norms[j]:violations+=1
    assert violations==0
    report=dict(passed=True,n=841,dimension=12,pairs_checked=841*840//2,
                interpretation='Exact printed decimal rays, normalized mathematically',
                certified_cosine_upper_numerator=p,certified_cosine_upper_denominator=q,
                margin_below_half='2^-24',source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                elapsed_seconds=time.monotonic()-start)
    (ROOT/'decimal841_certificate.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
