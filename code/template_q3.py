"""The template of Section 7 for q = 3, n = 3, r = 1: B_i = complement of the line S_i^perp in F_3^2."""
import itertools
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_srg import check_srg

P = list(itertools.product(range(3), repeat=2))
idx = {p: i for i, p in enumerate(P)}
lines = [frozenset(((a * d[0]) % 3, (a * d[1]) % 3) for a in range(3)) for d in [(1, 0), (0, 1), (1, 1), (1, 2)]]


def perp(L):
    return frozenset(v for v in P if all((v[0] * b[0] + v[1] * b[1]) % 3 == 0 for b in L))


def sub(x, y):
    return ((x[0] - y[0]) % 3, (x[1] - y[1]) % 3)


S = lines[:3]                                   # supports of the three pairs
A = np.zeros((27, 27), dtype=np.uint8)
for i in range(3):
    C = (perp(S[i - 1]) | perp(S[i])) - {(0, 0)}
    B = frozenset(P) - perp(S[i])               # H 1_B = 9 delta - 3 [beta in S_i]: z = -1 on S_i^*
    j = (i + 1) % 3
    for x in P:
        for y in P:
            if sub(x, y) in C:
                A[9 * i + idx[x], 9 * i + idx[y]] = 1
            if sub(y, x) in B:
                A[9 * i + idx[x], 9 * j + idx[y]] = A[9 * j + idx[y], 9 * i + idx[x]] = 1
check_srg(A, (27, 16, 10, 8))
print('q=3 template: srg(27,16,10,8) verified (the complement of the Schlaefli graph has parameters (27,10,1,5); '
      'srg(27,16,10,8) is unique, the Schlaefli graph)')
