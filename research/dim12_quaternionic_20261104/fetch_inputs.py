"""Restore public numerical inputs, checking the archived SHA256 hashes."""
from pathlib import Path
from urllib.request import urlretrieve
import argparse,hashlib,json

ROOT=Path(__file__).resolve().parent
URLS={
 'canonical840.txt':'https://raw.githubusercontent.com/k-nic/quaternionic/main/canonical840.txt',
 'gram841_coordinates.txt':'https://raw.githubusercontent.com/k-nic/quaternionic/main/gram841_coordinates.txt',
 'alignment_U.txt':'https://raw.githubusercontent.com/k-nic/quaternionic/main/canonical840_alignment_U.txt',
 'alignment_permutation.txt':'https://raw.githubusercontent.com/k-nic/quaternionic/main/canonical840_alignment_permutation.txt',
 'alignment_top24.txt':'https://raw.githubusercontent.com/k-nic/quaternionic/main/canonical840_alignment_top24.txt',
 'families.txt':'https://raw.githubusercontent.com/k-nic/quaternionic/main/canonical840_families.txt',
 'Gram842.txt':'https://raw.githubusercontent.com/k-nic/841_in_12D/main/beyond841/Gram842.txt',
 'Gram844.txt':'https://raw.githubusercontent.com/k-nic/841_in_12D/main/beyond841/Gram844.txt',
}


def main():
    p=argparse.ArgumentParser();p.add_argument('--include-near-misses',action='store_true');a=p.parse_args()
    pins=json.loads((ROOT/'inputs.json').read_text())['files'];data=ROOT/'data';data.mkdir(exist_ok=True)
    for name,url in URLS.items():
        if name.startswith('Gram') and not a.include_near_misses:continue
        path=data/name
        if not path.exists():
            temporary=path.with_suffix('.download')
            urlretrieve(url,temporary)
            if hashlib.sha256(temporary.read_bytes()).hexdigest()!=pins[name]['sha256']:
                raise RuntimeError(f'Public input changed: {name}; saved pins are not overwritten.')
            temporary.replace(path)
        assert hashlib.sha256(path.read_bytes()).hexdigest()==pins[name]['sha256'],name
        print('verified',name)


if __name__=='__main__':main()
