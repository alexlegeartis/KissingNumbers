"""Registry wrapper for the exact 67D signed support replay."""
from pathlib import Path
import importlib.util

def results():
    p=Path(__file__).resolve().parent/'verify.py'
    s=importlib.util.spec_from_file_location('signed67_replay',p)
    v=importlib.util.module_from_spec(s)
    s.loader.exec_module(v)
    return {67:v.main()['total']}

if __name__=='__main__':
    import json
    print(json.dumps(results(),sort_keys=True))
