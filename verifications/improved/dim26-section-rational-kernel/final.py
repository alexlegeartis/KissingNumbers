"""Registry entry for the independently replayable 26D section witness."""
import json
from pathlib import Path


def results():
    proof = json.loads(Path(__file__).with_name("verification.json").read_text())
    assert proof["status"] == "PASS" and proof["full_exterior_pair_check"]
    assert proof["fixed_exterior"] + proof["kernel"] == proof["cardinality"] == 201203
    return [{"dimension": proof["dimension"], "total": proof["cardinality"]}]
