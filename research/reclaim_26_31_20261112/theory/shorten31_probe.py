"""Reproduce the fixed-parent shortening and dual-weight audit."""
import argparse
import collections
import json
from pathlib import Path
import numpy as np
from source_utils import REPO

def words(G):
    masks=[sum(int(v)<<j for j,v in enumerate(row))for row in G]
    result=[0]
    for mask in masks:result += [x^mask for x in result]
    return result

def rref(G):
    G=G.copy();pivots=[];row=0
    for j in range(G.shape[1]):
        pivot=next((k for k in range(row,len(G))if G[k,j]),None)
        if pivot is None:continue
        G[[row,pivot]]=G[[pivot,row]]
        for k in range(len(G)):
            if k!=row and G[k,j]:G[k]^=G[row]
        pivots.append(j);row+=1
    return G[:row],pivots

def dual(G):
    R,pivots=rref(G);free=[j for j in range(G.shape[1])if j not in pivots]
    H=np.zeros((len(free),G.shape[1]),dtype=np.uint8)
    for i,j in enumerate(free):
        H[i,j]=1
        for k,q in enumerate(pivots):H[i,q]=R[k,j]
    assert not np.any((H@G.T)%2)
    return H

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',type=Path,default=REPO/'verifications/improved/dim32-35-signed-layers/source/public')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    G=np.array(json.loads((args.input_dir/'dense_code.json').read_text())['generator_matrix'],dtype=np.uint8)
    C=words(G);assert len(C)==131072 and min(x.bit_count()for x in C if x)==8
    supports=[int(x,16)for x in (args.input_dir/'a32.8.8.1713H').read_text().splitlines()if x and not x.startswith('$')]
    assert len(supports)==1713
    answer=[]
    for coordinate in range(32):
        shortened=G.copy();row=int(np.flatnonzero(shortened[:,coordinate])[0])
        for j in range(len(shortened)):
            if j!=row and shortened[j,coordinate]:shortened[j]^=shortened[row]
        shortened=np.delete(np.delete(shortened,row,axis=0),coordinate,axis=1)
        weights=dict(collections.Counter(x.bit_count()for x in words(dual(shortened))))
        remaining=[(x&((1<<coordinate)-1))|((x>>(coordinate+1))<<coordinate)for x in supports if not(x>>coordinate)&1]
        ranks=collections.Counter(len(rref(shortened[:,[j for j in range(31)if(x>>j)&1]])[1])for x in remaining)
        answer.append({'deleted':coordinate,'support_count':len(remaining),'proj_ranks':dict(ranks),'dual_weights':weights,'blind_count':65536+128*len(remaining)+1860})
    assert max(x['blind_count']for x in answer)==236996
    text=json.dumps(answer,indent=2)+'\n'
    if args.output:args.output.write_text(text)
    else:print(text,end='')

if __name__=='__main__':main()
