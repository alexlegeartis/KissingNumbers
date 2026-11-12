"""Prepare regenerable input files in an explicit scratch output directory."""
import argparse
import json
import shutil
from pathlib import Path
import numpy as np
from source_utils import ROOT,REPO,add_source_arguments,load_integer_source
from verify_small_certificates import golay_words

def main():
    parser=argparse.ArgumentParser(description=__doc__);add_source_arguments(parser)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--support',action='store_true',help='Prepare the complete signed-octad support pilot.')
    parser.add_argument('--layer',action='store_true',help='Prepare the 3+5 layer permutation annealer; requires the pinned 31D incumbent.')
    args=parser.parse_args();assert args.support or args.layer
    out=args.output_dir;out.mkdir(parents=True,exist_ok=True)
    if args.support:
        parent=REPO/'verifications/improved/dim32-35-signed-layers/source/public'
        shutil.copyfile(parent/'a32.8.8.1713H',out/'a32.8.8.1713H')
        shutil.copyfile(parent/'dense_code.json',out/'dense_code.json')
        (out/'golay_octads.hex').write_text(''.join(f'{w:x}\n'for w in golay_words()if w.bit_count()==8))
    if args.layer:
        _,_,X,_=load_integer_source(31,args)
        profiles,weights=np.unique(np.abs(X),axis=0,return_counts=True)
        moving=np.any(profiles[:,24:]!=0,axis=1)
        assert np.all(np.sort(profiles[~moving,:24],axis=1)[:,-3:].sum(1)<=16)
        profiles=profiles[moving];weights=weights[moving]
        assert np.max(np.count_nonzero(profiles[:,24:],axis=1))<=4
        with (out/'absolute_profiles31.txt').open('w')as stream:
            stream.write(f'{len(profiles)}\n')
            for row,weight in zip(profiles,weights):stream.write(' '.join(map(str,[int(weight),*map(int,row)]))+'\n')
        design=json.loads((ROOT/'layer_3plus5_design.json').read_text())
        with (out/'layer_3plus5_design.txt').open('w')as stream:
            for edge,triple in zip(design['tail_omitted_edges'],design['core_triples']):stream.write(' '.join(map(str,[edge[0]-1,edge[1]-1,*triple]))+'\n')
    print(json.dumps({'output_dir':str(out),'support_inputs':args.support,'layer_inputs':args.layer}))

if __name__=='__main__':main()
