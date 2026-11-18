"""Classification of all seeds on F_16^2 (Remark 3.7 of the paper).

Claim: every seed is z * chi_c for a C_5-invariant seed z and a character chi_c.

Reduction. The semilinear group GammaL(2,16) permutes seeds (Lemma 3.3 for GL(2,16); the Frobenius a -> a^2 applied
to both coordinates preserves Tr(ax + by)), and it maps C_5-invariant seeds to C_5-invariant seeds and character
twists to character twists. PGammaL(2,16) is 3-transitive on the 17 spread lines, so its orbits on four-sets of
spread lines are the cross-ratio classes of seed_search.py: 2040 four-sets of class (2,3,8,9,14,15) and 340 subline
four-sets of class (6,7). So it suffices to treat one four-set of each class:
    python newseeds.py 0,1,2,3 4      (non-subline class; no C_5-invariant seeds, so nothing is blocked)
    python newseeds.py 0,1,6,7 4      (the subline {0, 1, w, w^2}, w = t^2 + t)

Method. For the given four slopes, block every known seed (all C_5-invariant seeds on these slopes times all 256
characters) and solve the 0/1 model with CP-SAT: x_b in {0,1} gives z(b) = 2x_b - 1 on the 60 support points, and
y_v in {0,1} encodes Hz(v) = -6 + 16 y_v. INFEASIBLE means the known seeds are all seeds on this four-set.
python newseeds.py SLOPES WORKERS [TIMELIMIT]
"""
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
from ortools.sat.python import cp_model

sys.path.insert(0, str(Path(__file__).resolve().parent))
from srg768 import CHI, coset, is_seed, line   # noqa: E402
from seed_search import cross_ratio_class      # noqa: E402

slopes = [int(s) for s in sys.argv[1].split(',')]
workers = int(sys.argv[2])
print('slopes', slopes, 'cross-ratio class', cross_ratio_class(slopes), flush=True)
pts = [p for m in slopes for p in line(m)]
known = set()
for signs in itertools.product((1, -1), repeat=12):
    z = np.zeros(256, dtype=np.int64)
    for k, m in enumerate(slopes):
        for a in range(1, 16):
            z[line(m)[a - 1]] = signs[3 * k + coset(a)]
    if not is_seed(z):
        continue
    for c in range(256):
        known.add(tuple(int(z[b] * CHI[b, c]) for b in pts))
print('known seeds (C5-invariant times characters):', len(known), flush=True)
m = cp_model.CpModel()
x = [m.NewBoolVar(f'x{i}') for i in range(len(pts))]
y = [m.NewBoolVar(f'y{v}') for v in range(256)]
for v in range(256):
    coeffs = [int(CHI[b, v]) for b in pts]
    m.Add(sum(2 * c * xi for c, xi in zip(coeffs, x)) - 16 * y[v] == sum(coeffs) - 6)
for z in known:   # forbid exactly this sign vector: at least one literal differs
    m.AddBoolOr([x[i].Not() if zi == 1 else x[i] for i, zi in enumerate(z)])
solver = cp_model.CpSolver()
solver.parameters.num_workers = workers
solver.parameters.max_time_in_seconds = float(sys.argv[3]) if len(sys.argv) > 3 else 86400.0
t = time.time()
st = solver.Solve(m)
print('status', solver.StatusName(st), 'time', round(time.time() - t, 1), flush=True)
if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    z = [2 * solver.Value(xi) - 1 for xi in x]
    zz = np.zeros(256, dtype=np.int64)
    zz[pts] = z
    print('NEW SEED', is_seed(zz), json.dumps(z), flush=True)
