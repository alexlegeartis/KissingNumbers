from pathlib import Path
import numpy as np
p=Path(__file__).parent;z=np.load(p/'coordinates.npz');X=z['switched_twice'].astype(np.int32);assert np.all(X%18==0);(X//18).astype(np.int16).tofile(p/'search17_twice_source.bin')
