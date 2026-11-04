"""Full exact geometry and independent direct-integer changed-pair replay."""
from pathlib import Path
from itertools import product
import argparse,hashlib,json
import numpy as np
from verify_geometry import verify,read_points,require,canonical_hash
ROOT=Path(__file__).resolve().parent

def manifest():
    path=ROOT/'data/SHA256SUMS.txt'
    expected={name:sha for sha,name in (line.split('  ',1) for line in path.read_text().splitlines() if line)}
    actual={str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and
            '__pycache__' not in p.parts and p!=path and p.suffix!='.pyc'}
    require(set(expected)==actual,'Complete package manifest')
    for name,sha in expected.items():require(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha,'Hash: '+name)
    return len(expected)

def sign_cube31():
    rows=[]
    for a,b,c in product((-1,1),repeat=3):
        v=np.zeros(31,np.int16);v[[0,5,6,10]]=4*a;v[[8,20]]=4*b;v[24]=4*c;v[25]=-4*c
        rows.append(v)
    return np.array(rows,np.int16)

def changed_pairs(record):
    source=ROOT/record['baseline'];target=ROOT/record['certificate']
    require(hashlib.sha256(source.read_bytes()).hexdigest()==record['baseline_sha256'],'Baseline hash')
    require(hashlib.sha256(target.read_bytes()).hexdigest()==record['certificate_sha256'],'Certificate hash')
    old=read_points(source);final=read_points(target)
    delta=np.load(ROOT/record['delta'],allow_pickle=False);removed=delta['removed'];added=delta['added']
    for values in (removed,added):
        require(values.dtype.kind in 'iu' and values.ndim==2 and values.shape[1]==record['dimension'],'Delta shape')
        require(np.all(abs(values.astype(np.int64))<=16),'Delta coefficient bound')
        require(len(np.unique(values,axis=0))==len(values),'Distinct delta points')
        require(np.all(np.sum(values.astype(np.int64)**2,axis=1)==128),'Delta norms')
    oldset={row.tobytes() for row in old};finalset={row.tobytes() for row in final}
    rset={row.astype(np.int16).tobytes() for row in removed};aset={row.astype(np.int16).tobytes() for row in added}
    require(rset<=oldset and not aset&oldset,'Removed membership and new additions')
    require((oldset-rset)|aset==finalset and len(finalset)==len(final),'Exact reconstruction and distinctness')
    require(len(removed)==record['removed'] and len(added)==record['added'],'Delta counts')
    require(len(old)==record['previous_repository'] and len(final)==record['total'],'Recorded counts')
    require(len(added)-len(removed)==record['gain'],'Gain count')
    if record['dimension']==31:
        require(not len(removed) and aset=={row.tobytes() for row in sign_cube31()},'Exact eight-point sign cube')
        first=np.zeros(31,np.int16);first[[8,20]]=8
        second=np.zeros(31,np.int16);second[[24,25]]=[8,-8]
        require(all(row.tobytes() in oldset for row in (first,-first,second,-second)),'Four retained roots of the A3 completion')
    integers=final.astype(np.int64);lookup={row.tobytes():i for i,row in enumerate(final)};maximum=-128
    for row in added.astype(np.int64):
        dots=np.zeros(len(final),np.int64)
        for j in np.flatnonzero(row):dots+=integers[:,j]*int(row[j])
        dots[lookup[row.astype(np.int16).tobytes()]]=-128
        require(int(dots.max())<=64,'Direct integer inserted-point angle')
        maximum=max(maximum,int(dots.max()))
    return {'status':'PASS','dimension':record['dimension'],'points':len(final),'removed':len(removed),
            'added':len(added),'gain':record['gain'],'direct_integer_pairs':len(added)*(len(final)-1),
            'maximum_inserted_product':maximum,'canonical_point_sha256':canonical_hash(final)}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write-reports',action='store_true');args=parser.parse_args()
    if not args.write_reports:print(json.dumps({'manifest_files':manifest()}),flush=True)
    for record in json.loads((ROOT/'data/provenance.json').read_text())['records']:
        geometry=verify(ROOT/record['certificate']);delta=changed_pairs(record)
        require(geometry['points']==record['total'],'Full count')
        if args.write_reports:
            (ROOT/record['geometry_report']).write_text(json.dumps(geometry,indent=2)+'\n')
            (ROOT/record['delta_report']).write_text(json.dumps(delta,indent=2)+'\n')
        else:
            expected=json.loads((ROOT/record['geometry_report']).read_text())
            require({k:v for k,v in geometry.items() if k!='seconds'}=={k:v for k,v in expected.items() if k!='seconds'},'Fresh complete geometry matches report')
            require(delta==json.loads((ROOT/record['delta_report']).read_text()),'Fresh direct-integer proof matches report')
        print(json.dumps({'status':'PASS','geometry':geometry,'changed_pairs':delta}),flush=True)
if __name__=='__main__':main()
