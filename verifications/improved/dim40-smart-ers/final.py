"""Counting driver for the Smart / Edel–Rains–Sloane 40D construction."""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

def results():
    c = json.loads((HERE / "certificate.json").read_text())
    n, w, m = c["dimension"], c["support_weight"], c["supports"]
    assert (n, w, m, c["dense_words"]) == (40, 8, 4267, 589824)
    assert (c["support_minimum_distance"], c["dense_minimum_distance"]) == (8, 10)
    total = c["dense_words"] + m * (1 << (w - 1)) + 2 * n * (n - 1)
    assert total == c["total"] and total - c["previous_repository_bound"] == c["gain"]
    return [{"dimension": n, "total": total}]

if __name__ == "__main__":
    print(json.dumps(results()[0]))
