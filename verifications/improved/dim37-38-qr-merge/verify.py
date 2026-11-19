"""Run both complete exact layered certificates."""
import json
import verify37,verify38
if __name__=='__main__':
    print(json.dumps(verify37.replay(),sort_keys=True))
    verify38.main()
