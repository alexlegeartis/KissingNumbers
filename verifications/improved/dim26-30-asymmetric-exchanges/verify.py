"""Complete exact replay of the supplied, potentially non-antipodal codes."""
import argparse,hashlib,json
from pathlib import Path
from verify_core import verify
from verify_pool import verify as verify_pool
ROOT=Path(__file__).resolve().parent


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('dimensions',type=int,nargs='*');args=parser.parse_args()
    records=json.loads((ROOT/'data/provenance.json').read_text())['records'];available={r['dimension'] for r in records}
    assert set(args.dimensions)<=available
    for record in records:
        if args.dimensions and record['dimension'] not in args.dimensions:continue
        path=ROOT/record['certificate'];assert hashlib.sha256(path.read_bytes()).hexdigest()==record['certificate_sha256']
        result=verify(path);assert result['points']==record['total'] and result['dimension']==record['dimension']
        print(json.dumps(result,indent=2),flush=True)
        if record['dimension'] in (28,30):
            print(json.dumps(verify_pool(ROOT/f"data/finite_clique_cover{record['dimension']}.npz"),indent=2),flush=True)
