# Finer exact recurrence for the Γ72 rooted latitudes

The vertex sets, degree counts and per-vertex triangle caps are unchanged
from the [rooted exchange proof](../dim74-85-gamma72-rooted-exchange/).
We evaluate the same certified lower bound with up to 32,768 smaller
blocks instead of 512. The existing block inequality is valid for every
positive block length. Its conservatism is smaller on shorter blocks:
it bounds the recurrence coefficients over a shorter interval. This is
an arithmetic improvement, not a new tail geometry or new triangle cap.

| Dimension | Preceding | Verified lower bound | Gain |
|---:|---:|---:|---:|
| 74 | 6,218,213,672 | **6,218,213,857** | +185 |
| 75 | 6,218,331,258 | **6,218,332,102** | +844 |
| 76 | 6,218,636,373 | **6,218,639,446** | +3,073 |
| 77 | 6,219,014,956 | **6,219,020,876** | +5,920 |
| 78 | 6,219,721,604 | **6,219,734,842** | +13,238 |
| 79 | 6,220,599,336 | **6,220,622,144** | +22,808 |
| 80 | 6,221,725,998 | **6,221,762,278** | +36,280 |
| 81 | 6,222,202,604 | **6,222,229,444** | +26,840 |
| 82 | 6,224,435,960 | **6,224,479,692** | +43,732 |
| 83 | 6,227,388,086 | **6,227,455,178** | +67,092 |
| 84 | 6,231,161,058 | **6,231,259,216** | +98,158 |
| 85 | 6,235,855,650 | **6,235,993,850** | +138,200 |

The verifier pins and replays the parent degree and triangle calculations,
then independently checks the old 512-block result. Inside each new block
it performs one exact integer division and rounds the scaled recurrence
down. At the end the exact rational graph lower bound is rounded up to
an integer because the independence number is integral. Dimensions 74–76
use all body orientations, with exact degree multiplicities. Dimensions
77–85 give two antipodal output points per selected projective vertex.

Run `python verifications/improved/dim74-85-gamma72-fine-recurrence/verify.py`
with Python, NumPy, SymPy and assertions enabled. No external download or
numerical optimizer is used. The graph recurrence is due to
[Boßecker and Rautenbach](https://doi.org/10.1016/j.disc.2010.05.026).
