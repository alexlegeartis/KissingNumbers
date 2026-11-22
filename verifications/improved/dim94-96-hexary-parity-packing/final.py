"""Lower bound registry for the complementary six-row parity append."""
import json
from pathlib import Path
def main():
    return {r['dimension']:r['lower_bound'] for r in json.loads(Path(__file__).with_name('claims.json').read_text())}
if __name__=='__main__':
    print(main())
