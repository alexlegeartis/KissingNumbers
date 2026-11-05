from pathlib import Path
import numpy as np
z=np.load(Path(__file__).with_name('source.npz'))
np.c_[z['rational'],z['radical']].astype(np.int16).tofile('source.bin')
print('Wrote source.bin in current directory.')
