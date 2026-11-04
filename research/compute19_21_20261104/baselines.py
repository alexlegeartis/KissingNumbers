from pathlib import Path
import numpy as np

DATA=Path(__file__).resolve().parent/'data'


def baseline(dimension):
    z=np.load(DATA/f'baseline{dimension}.npz')
    if dimension==19:
        return np.concatenate([z['base']/np.sqrt(8.),z['axis']/np.sqrt(4889129.),
                               (z['diagR']+np.sqrt(6.)*z['diagS'])/np.sqrt(162.)])
    if dimension==21:
        return np.concatenate([z['base']/np.sqrt(8.),z['auxiliary']/np.sqrt(164.)])
    raise ValueError(dimension)
