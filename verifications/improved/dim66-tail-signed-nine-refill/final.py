"""Registry wrapper for the exact 66D signed refill."""
from pathlib import Path
import importlib.util

def results():
    p=Path(__file__).resolve().parent/'verify.py'
    s=importlib.util.spec_from_file_location('signed66_replay',p)
    v=importlib.util.module_from_spec(s)
    s.loader.exec_module(v)
    return {66:v.main()['total']}

if __name__=='__main__':
    import json
    print(json.dumps(results(),sort_keys=True))
