#!/usr/bin/env python3
"""Replay the pinned predecessor followed by the exact 93D long-shell verifier."""
import json
import shutil
import subprocess
import sys
from pathlib import Path

here = Path(__file__).resolve().parent
node = shutil.which("node") or shutil.which("nodejs")
if not node:
    raise SystemExit("Node.js 18 or newer is needed for the exact long-shell replay.")
predecessor = here.parent / "dim93-z8-split-d5-heads"
subprocess.run([sys.executable, str(predecessor / "verify.py")], check=True)
reported = json.loads(subprocess.check_output([node, str(here / "verify.mjs")], text=True))
claim = json.loads((here / "claims.json").read_text(encoding="utf-8"))
for key, observed in {"dimension": reported["dimension"], "previous": reported["old"],
                      "gain": reported["added"], "lower_bound": reported["lowerBound"]}.items():
    if claim[key] != observed:
        raise RuntimeError(f"Claim {key} does not match the exact witness")
print(json.dumps(reported, sort_keys=True))
