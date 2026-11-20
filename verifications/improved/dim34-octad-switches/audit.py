"""Exact 34D replay."""
import json
from replay import verify

if __name__ == '__main__':
    print(json.dumps(verify(expanded=True), sort_keys=True))
