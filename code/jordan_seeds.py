"""All generalized seeds on a four-set of type (c) (Proposition 9.2 of the paper).

Default: the support of the original seed (data/jordan_seed.json, F_2^8 with the dot product).
With 'normal': the four-set X, Y, D, J of Table 3 in V = F_4^4 with the pairing Tr(a_1 b_1 + ... + a_4 b_4)
(code/jordan_normal.py); the stored seed is then the seed of Table 3 (data/jordan_normal.json).

Iterated CP-SAT: find a seed with the given support, record its 256 character twists z * chi_c, exclude them, and
repeat until the model is INFEASIBLE. Prints the number of twist classes and checks that the stored seed is among them.
python jordan_seeds.py WORKERS [normal]
"""
import sys
import time
from pathlib import Path

import numpy as np
from ortools.sat.python import cp_model

sys.path.insert(0, str(Path(__file__).resolve().parent))
from jordan768 import CHAR, load        # noqa: E402

workers = int(sys.argv[1])
if 'normal' in sys.argv[2:]:
    import json
    from jordan_normal import CHN as CHAR    # noqa: F811
    z0 = np.array(json.loads((Path(__file__).resolve().parent.parent / 'data' / 'jordan_normal.json').read_text())['z'])
else:
    z0, _ = load()
pts = np.flatnonzero(z0).tolist()
m = cp_model.CpModel()
x = [m.NewBoolVar(f'x{i}') for i in range(len(pts))]
y = [m.NewBoolVar(f'y{v}') for v in range(256)]
for v in range(256):
    coeffs = [int(CHAR[v, b]) for b in pts]
    m.Add(sum(2 * c * xi for c, xi in zip(coeffs, x)) - 16 * y[v] == sum(coeffs) - 6)
classes, blocked = [], set()
start = time.time()
while True:
    solver = cp_model.CpSolver()
    solver.parameters.num_workers = workers
    st = solver.Solve(m)
    print('round', len(classes) + 1, solver.StatusName(st), 'time', round(time.time() - start, 1), flush=True)
    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        break
    z = np.zeros(256, dtype=np.int64)
    z[pts] = [2 * solver.Value(xi) - 1 for xi in x]
    assert set((CHAR @ z).tolist()) <= {10, -6}
    twists = {tuple(int(z[b] * CHAR[c, b]) for b in pts) for c in range(256)}
    assert len(twists) == 256 and not (twists & blocked)
    blocked |= twists
    classes.append(z)
    for t in twists:
        m.AddBoolOr([x[i].Not() if ti == 1 else x[i] for i, ti in enumerate(t)])
print('twist classes of seeds on this support:', len(classes), '; seeds:', len(blocked),
      '; stored seed among them:', tuple(int(z0[b]) for b in pts) in blocked, flush=True)
