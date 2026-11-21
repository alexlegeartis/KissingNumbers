"""Registry wrapper for the exact 61D and 62D replay."""
from pathlib import Path
import importlib.util

def results():
    p=Path(__file__).resolve().parent/'verify.py'
    s=importlib.util.spec_from_file_location('marked_replay',p)
    v=importlib.util.module_from_spec(s)
    s.loader.exec_module(v)
    return {d:x['total'] for d,x in v.main().items()}

if __name__=='__main__':
    import json
    print(json.dumps(results(),sort_keys=True))
