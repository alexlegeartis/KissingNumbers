from pathlib import Path
import numpy as np
z=np.load(Path(__file__).with_name('source.npz'));z['X'].astype(np.int16).tofile('source.bin')
print('Wrote source.bin in current directory.')
