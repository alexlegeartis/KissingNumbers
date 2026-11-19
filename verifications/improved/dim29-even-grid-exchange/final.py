import importlib.util
import json
from pathlib import Path

def results():
    path = Path(__file__).with_name("verify.py")
    spec = importlib.util.spec_from_file_location("verify_29d_even_grid_exchange", path)
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    _, checked = verifier.replay()
    return [{"dimension": checked["dimension"], "total": checked["total"]}]

if __name__ == "__main__":
    for row in results():
        print(json.dumps(row))
