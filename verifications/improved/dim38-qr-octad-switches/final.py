"""Return the 38D row after replaying the exact witness."""
import json
import importlib.util
from pathlib import Path


def replay():
    spec = importlib.util.spec_from_file_location("octad_switch_proof", Path(__file__).with_name("verify.py"))
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    return verifier.replay()


def results():
    result = replay()
    return [{"dimension": result["dimension"], "total": result["points"]}]


if __name__ == "__main__":
    for row in results():
        print(json.dumps(row))
