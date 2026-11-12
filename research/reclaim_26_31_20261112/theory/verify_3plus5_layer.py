"""Read-only exact verification of the 2,688-point layer and optional cost proof."""
import argparse
import itertools
import json
import numpy as np
from source_utils import ROOT,add_source_arguments,load_integer_source

def verify_geometry():
    design=json.loads((ROOT/'layer_3plus5_design.json').read_text())
    certificate=json.loads((ROOT/'layer_3plus5_verified.json').read_text())
    pc=certificate['core_permutation'];pt=certificate['tail_permutation']
    assert sorted(pc)==list(range(24)) and sorted(pt)==list(range(7))
    edges=design['tail_omitted_edges'];triples=design['core_triples']
    assert len(edges)==21 and sorted(map(tuple,edges))==list(itertools.combinations(range(1,8),2))
    supports=[]
    for edge,triple in zip(edges,triples):
        assert len(set(triple))==3 and all(0<=j<24 for j in triple)
        tail=sorted(set(range(7))-{pt[edge[0]-1],pt[edge[1]-1]})
        support=sorted([pc[j]for j in triple]+[24+j for j in tail])
        assert len(set(support))==8 and sum(j>=24 for j in support)==5
        supports.append(support)
    assert supports==certificate['supports']
    assert len(set(map(tuple,supports)))==21
    maximum=max(len(set(a)&set(b))for a,b in itertools.combinations(supports,2))
    assert maximum<=4
    parity=np.array([p for p in itertools.product((-1,1),repeat=8)if np.prod(p)==1],dtype=np.int16)
    assert len(parity)==128 and len(np.unique(parity,axis=0))==128
    minimum=int(min(np.count_nonzero(a!=b)for a,b in itertools.combinations(parity,2)))
    assert minimum==2
    Y=np.zeros((2688,31),dtype=np.int16)
    for i,support in enumerate(supports):Y[128*i:128*(i+1),support]=4*parity
    saved=np.load(ROOT/'layer_3plus5_points.npz',allow_pickle=False)
    assert np.array_equal(Y,saved['points'])
    assert len(np.unique(Y,axis=0))==2688
    assert np.all(np.sum(Y.astype(np.int64)**2,axis=1)==128)
    # Across supports, the product is <=16*intersection<=64.
    # Within a support it is <=16*(8-2*minimum_distance)=64.
    return supports,certificate,{'dimension':31,'new_layer_points':2688,'squared_norm':128,'distinct_pair_product_at_most':64,'maximum_support_intersection':maximum,'minimum_sign_distance':minimum,'exact_coordinate_reconstruction':True,'source_cost_checked':False}

def verify_cost(supports,certificate,args):
    filename,record,X,_=load_integer_source(31,args)
    assert np.max(np.count_nonzero(X[:,24:],axis=1))<=4
    # Every support has five tail positions; at least one is zero in every
    # incumbent vector. Its sign can enforce even parity at no product cost.
    deleted=np.flatnonzero(np.any([np.abs(X[:,s]).sum(1)>16 for s in supports],axis=0))
    saved=np.load(ROOT/'layer_3plus5_points.npz',allow_pickle=False)
    assert np.array_equal(deleted,saved['source_removals'])
    assert len(deleted)==certificate['deletion_count']==9649
    profiles,weights=np.unique(np.abs(X),axis=0,return_counts=True)
    B=np.array([profiles[:,s].sum(1)>16 for s in supports])
    assert int(weights[B.any(0)].sum())==len(deleted)
    relevant=np.flatnonzero(B.any(0));B=B[:,relevant];weights=weights[relevant]
    source=0;first_profile=22;sink=first_profile+len(relevant)
    balance=[0]*(sink+1);seen=set()
    for a,b,flow in certificate['positive_flow_edges']:
        assert (a,b)not in seen and flow>0;seen.add((a,b))
        if a==source and 1<=b<=21:capacity=128
        elif 1<=a<=21 and first_profile<=b<sink:
            assert B[a-1,b-first_profile];capacity=2689
        elif first_profile<=a<sink and b==sink:capacity=int(weights[a-first_profile])
        else:raise AssertionError('Invalid flow edge')
        assert flow<=capacity
        balance[a]+=flow;balance[b]-=flow
    assert balance[source]==2688 and balance[sink]==-2688
    assert not any(balance[1:sink])
    return {'source_filename':filename.name,'source_sha256':record['sha256'],'source_points':record['points'],'source_cost_checked':True,'shared_deletions':len(deleted),'full_layer_gain':2688-len(deleted),'maximum_gain_from_any_subset_of_these_21_whole_fibers':0,'integer_flow_checked':True,'scope':'The subset result concerns this one exhibited orientation and whole 128-sign fibers; partial sign fibers and other orientations are not excluded.'}

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--with-source',action='store_true',help='Also replay exact incumbent deletion costs and the integer max-flow certificate.')
    add_source_arguments(parser);args=parser.parse_args()
    supports,certificate,report=verify_geometry()
    if args.with_source:report.update(verify_cost(supports,certificate,args))
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
