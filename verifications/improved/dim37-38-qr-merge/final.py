"""Count the projected dense layers, octad families, and roots."""
from pathlib import Path
import json
import importlib.util
P=Path(__file__).resolve().parent
S=importlib.util.spec_from_file_location("qr_projective_supports",P/"supports37.py")
M=importlib.util.module_from_spec(S);S.loader.exec_module(M)
def results():
    res=[]
    for dim in (37,38):
        a=json.loads((P/f'dense_parent{dim+1}.json').read_text())
        n=(1<<len(a['generators_hex']))*len(a['representatives_hex'])
        if dim==38:n=3*n//4
        s=len(M.supports()) if dim==37 else len((P/'supports_3088.txt').read_text().splitlines())
        total=n+128*s+2*dim*(dim-1)
        assert total=={37:533864,38:643836}[dim]
        res.append({'dimension':dim,'total':total})
    return res
if __name__=='__main__':
    for row in results():
        print(json.dumps(row))
