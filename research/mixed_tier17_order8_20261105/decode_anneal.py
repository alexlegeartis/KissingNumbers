import numpy as np
p='research/mixed_tier17_order8_20261105/cache/';mask=np.fromfile(p+'anneal_potential.bin',np.uint8);assert len(mask)==15872 and np.all((mask==0)|(mask==1));chosen=np.flatnonzero(mask).astype(np.int32);assert len(chosen)>176
np.savez(p+'anneal_POTENTIAL.npz',chosen=chosen,kept=np.array([],np.int32));print('candidate',4322+8*len(chosen),'points: exact replay still required')
