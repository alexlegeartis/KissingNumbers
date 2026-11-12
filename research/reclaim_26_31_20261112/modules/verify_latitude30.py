"""Replay the finite norm8-core/norm24-tail obstruction in 30D.

Five explicit old blockers per candidate, together with the geometric
four-selected-candidates-per-old-point bound in README.md, exclude gain.
This is not a bound for any other candidate alphabet or kissing number.
"""
import argparse, hashlib, json
from itertools import combinations, product
from pathlib import Path
import numpy as np

I=np.asarray([0,2,4,5,6,10,11,23])

def catalogue():
    bodies=[]
    for a,b in combinations(range(8),2):
        for s,t in product((-2,2),repeat=2):
            row=np.zeros(8,np.int16);row[a]=s;row[b]=t;bodies.append(row)
    for signs in product((-1,1),repeat=8):
        if sum(s<0 for s in signs)%2==0:bodies.append(np.asarray(signs,np.int16))
    B=np.unique(np.asarray(bodies,np.int16),axis=0)
    T=np.asarray(list(product((-2,2),repeat=6)),np.int16)
    C=np.zeros((len(B)*len(T),30),np.int16);C[:,I]=np.repeat(B,len(T),axis=0);C[:,24:]=np.tile(T,(len(B),1))
    if len(B)!=240 or len(T)!=64 or len(np.unique(C,axis=0))!=15360:raise ValueError('Catalogue')
    if not np.all(np.sum(C.astype(np.int64)**2,axis=1)==32):raise ValueError('Candidate norms')
    return B,T,C

def main():
    p=argparse.ArgumentParser();p.add_argument('source',type=Path)
    p.add_argument('--certificate',type=Path,default=Path(__file__).parent/'data/latitude30_five_blockers.npz')
    a=p.parse_args();z=np.load(a.source,allow_pickle=False);X=z['points'].astype(np.int16)
    if int(z['squared_norm'])==128:
        if np.any(X%2):raise ValueError('Source scale')
        X//=2
    w=np.load(a.certificate,allow_pickle=False);expected=str(w['source_sha256'].item())
    if hashlib.sha256(a.source.read_bytes()).hexdigest()!=expected:raise ValueError('Source hash')
    B,T,C=catalogue();owners=w['owners']
    if owners.shape!=(15360,5) or np.any(owners<0) or np.any(owners>=len(X)):raise ValueError('Owner shape')
    if np.any(np.diff(np.sort(owners,axis=1),axis=1)==0):raise ValueError('Distinct five owners')
    dots=np.sum(C[:,None,:].astype(np.int64)*X[owners].astype(np.int64),axis=2)
    if np.any(dots<=16):raise ValueError('Forbidden products')
    old=X[np.unique(owners)]
    tails=old[:,24:]
    if np.any(np.count_nonzero(tails,axis=1)!=4) or np.any((tails!=0)&(abs(tails)!=2)):raise ValueError('Four-sign old tails')
    if np.max(B.astype(np.int64)@old[:,I].astype(np.int64).T)>8:raise ValueError('Old core product bound')
    if not np.all(np.sum(old.astype(np.int64)**2,axis=1)==32):raise ValueError('Old point norms')
    report={'status':'PASS','candidates':len(C),'explicit_forbidden_products':int(owners.size),
        'minimum_witness_product':int(dots.min()),'candidate_blocker_degree':5,
        'maximum_selected_conflicts_per_deleted_old_point':4,
        'conclusion':'For every selected subset, 5*added <= 4*deleted; no positive gain.',
        'scope':'Only the explicitly regenerated 15360-point norm8 E8-body/norm24 signed-tail family.'}
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
