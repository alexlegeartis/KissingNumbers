"""Bounded 38D pilot: full fourteen-coordinate sign heads on two unused roots.

At body fraction 1/4 these are automatically compatible with the retained grid
and Hamming layer, but cap and flat-axis deletions prove costly. This is a search
diagnostic, not an impossibility theorem for general latitude constructions.
"""
from pathlib import Path
import importlib.util
import json
import numpy as np

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]
PACKAGE=REPO/'verifications/improved/dim38-joint-latitudes'


def main():
    spec=importlib.util.spec_from_file_location('retained38_coordinates',PACKAGE/'refill/verify_coordinates.py')
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    A,_=old.reconstruct()
    source=PACKAGE/'source/source'
    labels=np.loadtxt(source/'class_labels.txt',dtype=np.int32)
    W=np.loadtxt(source/'dim14_directions.txt',delimiter=',',dtype=np.int16)
    triples=np.loadtxt(source/'triples.txt',dtype=np.int32)
    owner=np.empty(1932,np.int32)
    for i,row in enumerate(triples):owner[row]=i
    c=json.loads((PACKAGE/'refill/certificate.json').read_text())
    Q=np.array(c['rotation_numerator'],np.int16)
    axes=W[c['axis_keep']]@Q.T
    # First sign fixed positive: 8,192 antipodal head lines.
    signs=(1-2*((2*np.arange(8192,dtype=np.int32)[:,None]>>np.arange(14))&1)).astype(np.int16)
    roots=np.zeros((2,24),np.int16);roots[:,22]=4;roots[:,23]=[4,-4]
    contacts=np.array([np.bincount(labels[np.abs(A@u)==16],minlength=644)for u in roots])
    cost=np.empty((2,len(signs)),np.int32)
    axis_cost=np.empty(len(signs),np.int16)
    capmax=np.empty(len(signs),np.int16);axismax=capmax.copy()
    for start in range(0,len(signs),256):
        block=signs[start:start+256]
        capdot=np.abs(block@W.T);axisdot=np.abs(block@axes.T)
        capmax[start:start+len(block)]=capdot.max(axis=1)
        axismax[start:start+len(block)]=axisdot.max(axis=1)
        cost[:,start:start+len(block)]=2*contacts[:,owner]@(capdot==8).T
        axis_cost[start:start+len(block)]=(axisdot==8).sum(axis=1)
    report={'antipodal_head_lines_screened':len(signs),
            'uniform_cap_axis_survivors':int(np.sum((capmax<=6)&(axismax<=6))),
            'minimum_additional_cap_point_deletions_per_body_line':cost.min(axis=1).tolist(),
            'axis_deletion_histogram':{str(k):int(np.sum(axis_cost==k))for k in sorted(set(axis_cost))},
            'candidate_new_points_per_head_and_body_line':4,
            'scope':'Four independent-sign points on each tested head/body pair. Cap cost excludes its own body line. No joint-union optimality claim.',
            'body_roots':roots.tolist()}
    (HERE/'dense14_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
