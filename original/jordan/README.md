# Original scripts for the Jordan Cayley graph

These scripts first produced the Cayley graph of Section 9 of the paper (18 November 2026). They are copied unchanged
except for one line: in `local.py`, `OLD` now points at `..` (this repository's `original/`, which holds the same
`construct768.py`) instead of the original sibling directory `construction_routes_20261118`. Numpy is required.

| file | contents |
|---|---|
| `cayley_control.py` | searches the stored spreads for three Jordan-type four-sets cycled by a linear map M of order 3 with one fixed point, assembles the graph, checks all pairs, and writes the graph and its certificate |
| `assemble.py` | signed cross blocks on a binary partial spread (`Ports`), imported by `cayley_control.py` |
| `local.py` | field helpers and the type of a four-set of 4-spaces (semisimple of order 3, nonsemisimple of order 6, irreducible) |
| `crossratio.py`, `crossratio.json` | the twelve sign classes, up to translation, on each type of four-set; none on the irreducible type |
| `planes.py`, `planes.json` | spreads of order 16 from semifield multiplications, net replacement and duality; the graph uses `semifield_0` |
| `cayley_control.json` | the certificate written by `cayley_control.py` (identical to `../../data/jordan_seed.json`) |
| `group.g` | GAP script identifying the group from the basis images of R |

`verify_cayley.py` and `compare.py` from the original directory are not included: the first imports a verifier from
unrelated unpublished data, and the second compares against local archives. Their checks are repeated by
`../../code/jordan768.py`.

Rebuild and compare (all byte-identical, checked 18 November 2026):

```sh
python -B crossratio.py > /tmp/crossratio.json && cmp /tmp/crossratio.json crossratio.json
python -B planes.py > /tmp/planes.json && cmp /tmp/planes.json planes.json
python -B cayley_control.py /tmp/jordan.g6 > /tmp/seed.json
cmp /tmp/jordan.g6 ../../graphs/srg768_312_136_120_jordan_cayley.g6
cmp /tmp/seed.json cayley_control.json
```
