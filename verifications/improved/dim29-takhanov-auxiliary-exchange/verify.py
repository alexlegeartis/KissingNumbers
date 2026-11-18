"""Verify the pinned two-root and auxiliary exchange in dimension 29."""

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CHECKER = ROOT.parent / "dim30-contact-sum-exchanges" / "verify.py"


if __name__ == "__main__":
    spec = importlib.util.spec_from_file_location("integer_exchange_verifier", CHECKER)
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    verifier.main(ROOT)
