"""Replay the exact 370,416-point witness."""
import json
from replay import replay

if __name__ == '__main__':
    print(json.dumps(replay(expanded=True), sort_keys=True))
