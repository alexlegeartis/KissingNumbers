"""Reconstruct all coordinates from the frozen baseline and 4-block formula."""
import argparse,json,numpy as np
from pathlib import Path
from construction import ROOT,reconstruct,require
from verify_geometry import canonical_hash
p=argparse.ArgumentParser();p.add_argument('--output',type=Path);a=p.parse_args()
x=reconstruct();expected=np.load(ROOT/'data/coordinates31_251200.npz',allow_pickle=False)['points']
require(np.array_equal(x,expected),'Complete coordinate reconstruction')
if a.output:np.savez_compressed(a.output,points=x,squared_norm=128,inner_product_limit=64,coordinate_denominator=2)
print(json.dumps({'status':'PASS','dimension':31,'points':len(x),'canonical_point_sha256':canonical_hash(x)}))
