#!/usr/bin/env python3
"""Run the exact 60D Node verifier and record its pinned inputs for run_all.py."""
from pathlib import Path
import subprocess

here = Path(__file__).resolve().parent
root = here.parents[2]
for rel in ("verifications/improved/dim60-61-rm80-puncture/extra60.json.gz.b64", "verifications/improved/dim60-61-rm80-guard-transfer/transfer60.txt.gz.b64", "verifications/improved/dim62-67-rm2-cubic-cosets/witness.txt",):
    (root / rel).read_bytes()  # record cross-package inputs in the replay trace
subprocess.run(["node", str(here / "verify.mjs"), str(root)], check=True)
