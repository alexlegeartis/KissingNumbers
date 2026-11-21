"""Result adapter (exact replay) for the 86D ternary chord section."""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from verify import results as verify_results
def results():
    return [dict(dimension=int(d),total=item["total"],previous=item["previous"],gain=item["gain"])
            for d,item in verify_results().items()]
if __name__=="__main__":
    import json
    print(json.dumps(results(),indent=2))
