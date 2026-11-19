"""Return the 37D radical-layer count."""
import json
from pathlib import Path
P=Path(__file__).resolve().parent
def results():
    n=len((P/'six_supports.txt').read_text().split())
    return [{'dimension':37,'total':533864-144+64*n}]
if __name__=='__main__':
    for row in results():
        print(json.dumps(row))
