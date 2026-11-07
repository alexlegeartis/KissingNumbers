"""Demonstrate that deleted/replaced source points cannot simply be kept."""
from pathlib import Path
import json
import numpy as np

HERE=Path(__file__).resolve().parent
z=np.load(HERE/'verified_delta_coordinates.npz');G=z['gram_36']
a=z['new_equator']@G@z['old_base'].T
b=np.concatenate([z['new_base'],z['new_equator']])@G@z['old_equator'].T
assert a.max()>6 and b.max()>6
r={
    'restore_original_base_switches':{
        'detected_conflict':True,
        'new_equator_vs_original_base_max':int(a.max()),
        'first_max_indices':list(map(int,np.unravel_index(np.argmax(a),a.shape)))},
    'restore_removed_phase_points':{
        'detected_conflict':True,
        'changed_or_added_vs_removed_phase_max':int(b.max()),
        'first_max_indices':list(map(int,np.unravel_index(np.argmax(b),b.shape)))}}
(HERE/'negative_controls.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r,indent=2))
