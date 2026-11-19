"""Run both exact verification implementations."""
from pathlib import Path
import subprocess
import sys

if not __debug__:
    raise SystemExit('Run without -O: verification assertions must be enabled.')
root = Path(__file__).resolve().parent
for script in ['verify_shared.py', 'audit_coordinates.py']:
    subprocess.run([sys.executable, str(root / script)], cwd=root, check=True)
