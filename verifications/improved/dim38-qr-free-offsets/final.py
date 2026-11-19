"""Count the three dense families and two retained sparse layers."""
from pathlib import Path
P=Path(__file__).resolve().parent
def results():
    n=sum(int(s.split()[0]) for s in (P/'offset_witness.txt').read_text().splitlines())
    s=len((P.parent/'dim37-38-qr-merge'/'supports_3088.txt').read_text().splitlines())
    total=n*8192+s*128+2*38*37
    assert total==652028
    return [{'dimension':38,'total':total}]
if __name__=='__main__':
    import json
    for row in results():
        print(json.dumps(row))
