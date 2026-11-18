# Original construction scripts

These are the scripts that first produced the 768-vertex graphs (18 November 2026), copied unchanged except where
noted. The code in `../code/` was written separately from the mathematical description and does not import them.

| file | contents |
|---|---|
| `construct768.py` | builds the three-orbit (shear) graph, its complement and the Cayley graph; checks all pairs, counts four-cliques, verifies the Cayley connection set |
| `construct_405.py` | helper module imported by `construct768.py` (its `check` and `graph6` functions); copied from the 405 construction directory |
| `ports768.py`, `ports768.jsonl` | the seed search over the two semilinear orbits of four spread lines, and its output (12 signings on one orbit, none on the other) |
| `README_construction_routes.md` | the original write-up. Links were changed to point at `../graphs/` and `../code/verify_srg.py`; a link to unrelated search notes was removed |
| `jordan/` | the scripts that produced the Jordan-type Cayley graph of Section 9 (`../graphs/srg768_312_136_120_jordan_cayley.g6`); see `jordan/README.md` |

The graph files written by `construct768.py` were renamed in `../graphs/`:

| original name | name here |
|---|---|
| `srg768_312_136_120.g6` | `srg768_312_136_120_threeorbit.g6` |
| `srg768_312_136_120_cayley.g6` | `srg768_312_136_120_cayley.g6` |
| `srg768_455_262_280.g6` | `srg768_455_262_280_threeorbit_complement.g6` |

Rebuild (standard library only), then compare:

```sh
python -B construct768.py --output /tmp/rebuild768
cmp /tmp/rebuild768/srg768_312_136_120.g6 ../graphs/srg768_312_136_120_threeorbit.g6
cmp /tmp/rebuild768/srg768_312_136_120_cayley.g6 ../graphs/srg768_312_136_120_cayley.g6
cmp /tmp/rebuild768/srg768_455_262_280.g6 ../graphs/srg768_455_262_280_threeorbit_complement.g6
```

All three are byte-identical (checked 18 November 2026).
