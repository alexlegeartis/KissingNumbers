"""Resolve explicitly supplied, hash-pinned incumbent files."""
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[2]

def add_source_arguments(parser):
    parser.add_argument('--source-dir',action='append',type=Path,default=[],help='Directory containing pinned NPZ sources; repeat for multiple directories. A PackingStar repository root is also accepted.')
    parser.add_argument('--source',action='append',default=[],metavar='DIM=FILE',help='Explicit source-file override; contents must match the pinned hash.')

def resolve_source(dimension,args):
    record=json.loads((ROOT/'source_manifest.json').read_text())['sources'][str(dimension)]
    overrides={}
    for item in args.source:
        dim,filename=item.split('=',1);overrides[int(dim)]=Path(filename)
    candidates=[]
    if dimension in overrides:candidates=[overrides[dimension]]
    else:
        for base in args.source_dir:
            candidates += [base/record['filename'],base/'new_bounds_20261111'/record['filename'],base/'25D-31D/new_bounds_20261111'/record['filename']]
    for filename in candidates:
        if filename.is_file():
            actual=hashlib.sha256(filename.read_bytes()).hexdigest()
            if actual!=record['sha256']:raise ValueError(f'Wrong source hash: {filename}; expected {record["sha256"]}, got {actual}')
            return filename,record
    raise FileNotFoundError(f'Missing pinned source for dimension {dimension}: {record["filename"]}. Supply --source-dir DIRECTORY or --source {dimension}=FILE. See README.md.')

def load_integer_source(dimension,args):
    filename,record=resolve_source(dimension,args)
    data=np.load(filename,allow_pickle=False)
    if dimension==26:
        raw=data['rational_points'].astype(np.int64);den=int(data['rational_denominator'])
        if den<=0 or np.any((2*raw)%den):raise ValueError('Source cannot be expressed in the supplied integer scale')
        points=2*raw//den
        if len(points)+len(data['quadratic_A'])!=record['points']:raise ValueError('Wrong source population')
    else:
        points=data['points'].astype(np.int64);norm=int(data['squared_norm'])
        if norm==32:points*=2
        elif norm!=128:raise ValueError(f'Unexpected source norm: {norm}')
        if len(points)!=record['points']:raise ValueError('Wrong source population')
    if not np.all(np.sum(points*points,axis=1)==128):raise ValueError('Wrong source point norm')
    return filename,record,points,data

def radical_positive(a,b,d):
    """Exact sign of a+b*sqrt(d), for integer a,b and positive integer d."""
    a,b,d=int(a),int(b),int(d)
    if d<=0:raise ValueError('Radicand must be positive')
    if b==0:return a>0
    if a==0:return b>0
    if (a>0)==(b>0):return a>0
    difference=a*a-d*b*b
    return difference!=0 and ((difference>0)==(a>0))
