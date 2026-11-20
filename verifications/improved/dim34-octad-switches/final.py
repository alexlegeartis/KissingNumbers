"""Return the 34D row after exact reconstruction and verification."""
import importlib.util
import json
from pathlib import Path


def results():
    spec = importlib.util.spec_from_file_location("dim34_switch_replay", Path(__file__).with_name("replay.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return [{"dimension": 34, "total": module.verify()["points"]}]


if __name__ == "__main__":
    for row in results():
        print(json.dumps(row))
