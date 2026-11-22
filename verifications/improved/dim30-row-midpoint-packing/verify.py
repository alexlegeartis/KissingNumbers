"""Verify the row midpoint packing by exact exchange and all-pairs checks."""

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parent
previous_verifier = ROOT.parent / "dim30-contact-sum-exchanges" / "verify.py"
spec = importlib.util.spec_from_file_location("exchange_verifier", previous_verifier)
verifier = importlib.util.module_from_spec(spec)
spec.loader.exec_module(verifier)

if __name__ == "__main__":
    verifier.main(package=ROOT)
