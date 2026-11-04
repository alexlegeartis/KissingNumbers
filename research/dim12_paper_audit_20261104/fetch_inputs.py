"""Restore the public paper, code and numerical inputs with SHA256 checks."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import hashlib,json
from pathlib import Path
from urllib.request import urlretrieve
ROOT=Path(__file__).resolve().parent
URLS={
 'paper.pdf':'https://arxiv.org/pdf/2606.18984v2',
 'saved_48_systems.zip':'https://raw.githubusercontent.com/k-nic/841_in_12D/main/saved_48_systems.zip',
 'search_841_riesz.py':'https://raw.githubusercontent.com/k-nic/841_in_12D/main/search_841_riesz.py',
 'polish_841.py':'https://raw.githubusercontent.com/k-nic/841_in_12D/main/polish_841.py',
 'checking.ipynb':'https://raw.githubusercontent.com/k-nic/841_in_12D/main/checking.ipynb',
 'gram841_coordinates.txt':'https://raw.githubusercontent.com/k-nic/841_in_12D/main/gram841_coordinates.txt',
}


def main():
    pins=json.loads((ROOT/'audit.json').read_text())['source_files'];out=ROOT/'inputs';out.mkdir(exist_ok=True)
    for name,url in URLS.items():
        p=out/name
        if not p.exists():
            temp=p.with_suffix('.download');urlretrieve(url,temp)
            assert hashlib.sha256(temp.read_bytes()).hexdigest()==pins[name]['sha256'],name
            temp.replace(p)
        assert hashlib.sha256(p.read_bytes()).hexdigest()==pins[name]['sha256'],name
        print('verified',name)


if __name__=='__main__':main()
