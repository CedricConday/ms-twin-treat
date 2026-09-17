import json

import backtest.clinical as C
import bricks.readout as R

for scale in (1.5, 0.8, 3.0):
    R.MAX_RELAPSE = scale
    g = C.run_gate()
    print("=== MAX_RELAPSE =", scale, "| type:", type(g).__name__, flush=True)
    print(json.dumps(g, default=str)[:1200], flush=True)
