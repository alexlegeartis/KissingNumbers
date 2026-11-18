# Joint midpoint exchanges in dimensions 30 and 31

The full coordinate witnesses give **248,289 in 30D** and **251,524 in 31D**. The sources are the preceding repository certificates at revision [`fe5d9f7`](https://github.com/alexlegeartis/KissingNumbers/tree/fe5d9f7f63c7fc5dbe209e6bdc7693fbab58dbfc).

| Dimension | Source | Inserted | Removed | Net gain | Witness count | Witness unordered pairs | Current repository bound |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 30 | 248,212 | 545 | 468 | 77 | **248,289** | 30,823,589,616 | 248,351 |
| 31 | 251,420 | 274 | 170 | 104 | **251,524** | 31,632,035,526 | 251,524 |

The 248,289-point witness is the source for the [30D contact-sum continuation](../dim30-contact-sum-exchanges/).

The search used row midpoints between Leech vectors with a common octad projection. If the endpoints have squared norm 128 and product 64, their midpoint has squared norm 96. A tail of squared norm 32 restores the target squared norm. Reusing endpoints already omitted from the source makes these insertions cheaper. The joint exchanges also change fibre midpoints, auxiliary points and retained Leech rows. The 30D search continued with candidates obtained by summing points of inner product −64 and optimizing the resulting exchanges.

## Replay

Run from the repository clone with Python and NumPy.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verify.py
python final.py
```

`verify.py` checks the source, coordinate and exchange hashes, reconstructs each final set from its pinned source, and independently checks every pair involving an inserted point in signed int64. The shared [integer geometry checker](../dim28-takhanov-root-exchange/verify.py) then checks the distinctness, norms and **every unordered pair** in the final configurations. The checks rely on the actual stored coordinates. Each valid coordinate has magnitude at most 11, so all integer partial sums are smaller than `31 * 11**2 < 2**24` and binary32 matrix products are exact. Use `--direct` for the slower signed-int64 full replay, `--dimension 30` or `--dimension 31` for one record, and `--report PATH` for a fresh JSON receipt.

The `data/coordinates*.npz` files contain the complete integer witnesses. The `data/exchange*.npz` files hold the added and removed rows. `manifest.json` pins file identities and canonical coordinate hashes. Normalize each coordinate row by dividing by `sqrt(128)`.
