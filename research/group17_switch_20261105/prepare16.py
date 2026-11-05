from pathlib import Path
import numpy as np
p=Path(__file__).parent;z=np.load(p/'coordinates.npz');X=z['switched_twice'].astype(np.int32);X=X[X[:,16]==0,:16];assert X.shape==(4320,16) and np.all(X%108==0);X=(X//108).astype(np.int8);assert np.all((X.astype(np.int32)**2).sum(1)==8);X.tofile(p/'search16_switched_twice_baseline.bin');print('prepared 4320-point switched equator')
