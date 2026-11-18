"""Normal form of the Jordan-type Cayley graph (Section 8 of the paper, Table 3).

V = F_4^4 with F_4 = {0, 1, w, w^2}, encoded with two bits per coordinate (bit 2i: coefficient of 1, bit 2i+1:
coefficient of w in a_{i+1}). Pairing <a, b> = Tr_{F_4/F_2}(a_1 b_1 + ... + a_4 b_4), Tr(w) = Tr(w^2) = 1.
Support of the seed: X = {(u,0)}, Y = {(0,u)}, D = {(u,u)}, J = {(u, Lu)}, u in F_4^2, L(p,q) = (wp + q, wq).

python jordan_normal.py

1. Transports the seed and the map M of data/jordan_seed.json (F_2^8 with the dot product) to this normal form:
   P maps the four support subspaces to X, Y, D, J (all 24 orders, all solutions K of K L = L' K), M' = P M P^-1,
   z' = z o P^-1. Among the choices it takes the one where M'(x) = A x + B xbar has fewest nonzero entries in
   A, B (first in enumeration order), then the character twist z' chi_c with fewest minus signs (smallest c).
   Writes data/jordan_normal.json.
2. Builds the graph of Construction 4.1 from (z', M') with the pairing above, independently of the transport:
   R = adjoint inverse of M', B_i = {H z_i = 10}, C_i = punctured annihilators; checks it is
   srg(768,312,136,120) on all pairs, that tau is an automorphism, and that it is isomorphic to
   graphs/srg768_312_136_120_jordan_cayley.g6 (bliss canonical hash).
3. Prints Table 3 (signs) and the matrices A, B (called Q_1, Q_2 in the paper).
"""
import functools
import json
import sys
from itertools import combinations, permutations, product
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import jordan768 as J                              # noqa: E402
from verify_srg import check_srg, read_graph6      # noqa: E402

MUL = [[0, 0, 0, 0], [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2]]   # F_4: 0, 1, w, w^2 = 0, 1, 2, 3
CONJ = [0, 1, 3, 2]
TR = [0, 0, 1, 1]
NAMES = ['0', '1', 'w', 'w2']


def enc(a):
    return a[0] | a[1] << 2 | a[2] << 4 | a[3] << 6


def dec(x):
    return [(x >> (2 * i)) & 3 for i in range(4)]


def pair(x, y):
    s = 0
    for a, b in zip(dec(x), dec(y)):
        s ^= MUL[a][b]
    return TR[s]


DOTN = np.array([[pair(x, y) for y in range(256)] for x in range(256)], dtype=np.int64)
CHN = 1 - 2 * DOTN


def Lmap(p, q):
    return MUL[2][p] ^ q, MUL[2][q]


U = [(p, q) for q in range(4) for p in range(4) if (p, q) != (0, 0)]
SPACES = {
    'X': [enc((p, q, 0, 0)) for p, q in U],
    'Y': [enc((0, 0, p, q)) for p, q in U],
    'D': [enc((p, q, p, q)) for p, q in U],
    'J': [enc((p, q) + Lmap(p, q)) for p, q in U],
}


def AB_of(Mn):
    """M(x) = A x + B conj(x) over F_4: M(e_j) = A_j + B_j and M(w e_j) = w A_j + w^2 B_j."""
    A = [[0] * 4 for _ in range(4)]
    B = [[0] * 4 for _ in range(4)]
    for j in range(4):
        m1 = dec(int(Mn[enc([int(k == j) for k in range(4)])]))
        mw = dec(int(Mn[enc([2 * int(k == j) for k in range(4)])]))
        for i in range(4):
            b = mw[i] ^ MUL[2][m1[i]]
            A[i][j], B[i][j] = m1[i] ^ b, b
    return A, B


def apply_AB(A, B, x):
    v = dec(x)
    return enc([functools.reduce(lambda s, k: s ^ MUL[A[i][k]][v[k]] ^ MUL[B[i][k]][CONJ[v[k]]], range(4), 0)
                for i in range(4)])


def transport():
    z, M = J.load()
    sup = J.support_subspaces(z)
    Mm = np.array([J.apply(M, x) for x in range(256)])
    L0 = [enc(list(Lmap(*[(u >> 2 * i) & 3 for i in range(2)])) + [0, 0]) & 15 for u in (1, 2, 4, 8)]
    best = None
    for order in permutations(range(4)):
        A_, B_, C_, D_ = [sup[k] for k in order]
        a, b = J.basis_of(A_), J.basis_of(B_)
        co = J.solve_coords(a + b)

        def graph_map(S):
            m = {}
            for c in S:
                k = co[c]
                m[k & 15] = k >> 4
            return m
        phi, psi = graph_map(C_), graph_map(D_)
        phinv = {v: k for k, v in phi.items()}
        Lim = [phinv[psi[1 << j]] for j in range(4)]
        for Kim in product(range(1, 16), repeat=4):
            Kim = list(Kim)
            if len({J.apply(Kim, x) for x in range(16)}) < 16:
                continue
            if not all(J.apply(Kim, L0[j]) == J.apply(Lim, Kim[j]) for j in range(4)):
                continue
            # P^-1 (u, w) = K u + phi(K w): X -> A_, Y -> B_, D -> graph(phi) = C_, J -> graph(phi L) = D_
            Pinv = np.array([J.apply(a, J.apply(Kim, x & 15)) ^ J.apply(b, phi[J.apply(Kim, x >> 4)])
                             for x in range(256)])
            P = np.argsort(Pinv)
            Mn = P[Mm[Pinv]]
            A, B = AB_of(Mn)
            nz = sum(1 for r in A + B for e in r if e)
            if best is None or nz < best[0]:
                best = (nz, Mn, A, B, z[Pinv])
    nz, Mn, A, B, zn = best
    negs = [int(((zn * CHN[c]) < 0).sum()) for c in range(256)]
    c = int(np.argmin(negs))
    zt = zn * CHN[c]
    assert all(apply_AB(A, B, x) == Mn[x] for x in range(256))
    return {'A': A, 'B': B, 'z': [int(v) for v in zt], 'twist': c}


def build_normal(z, A, B):
    z = np.array(z)
    Mn = np.array([apply_AB(A, B, x) for x in range(256)])
    assert len(set(Mn.tolist())) == 256
    assert (Mn[Mn[Mn]] == np.arange(256)).all() and int((Mn == np.arange(256)).sum()) == 1
    assert set(np.flatnonzero(z).tolist()) == {x for S in SPACES.values() for x in S}
    zs = [z]
    for _ in range(2):
        nz = np.zeros(256, dtype=np.int64)
        nz[Mn] = zs[-1]                              # z_{i+1} = z_i o M^-1
        zs.append(nz)
    subs = [[frozenset([0] + [int(Mn[Mn[x]]) if i == 2 else int(Mn[x]) if i == 1 else x for x in S])
             for S in SPACES.values()] for i in range(3)]
    flat = [W for s in subs for W in s]
    assert all(len(X & Y) == 1 for X, Y in combinations(flat, 2)), 'the twelve subspaces must meet pairwise in 0'
    Bs = []
    for zi in zs:
        Hz = CHN @ zi
        assert set(Hz.tolist()) == {10, -6}
        Bs.append(np.flatnonzero(Hz == 10))
    Cs = []
    for i in range(3):
        C = set()
        for W in subs[i - 1] + subs[i]:
            C |= {x for x in range(1, 256) if not any(DOTN[x, w] for w in W)}
        assert len(C) == 120
        Cs.append(sorted(C))
    Adj = np.zeros((768, 768), dtype=np.uint8)
    x = np.arange(256)
    for i in range(3):
        for c in Cs[i]:
            Adj[256 * i + x, 256 * i + (x ^ c)] = 1
        j = (i + 1) % 3
        for b in Bs[i]:
            Adj[256 * i + x, 256 * j + (x ^ b)] = 1
            Adj[256 * j + (x ^ b), 256 * i + x] = 1
    # R = adjoint inverse of M: <R x, M b> = <x, b>
    R = np.array([next(y for y in range(256) if all(DOTN[y, Mn[1 << k]] == DOTN[xx, 1 << k] for k in range(8)))
                  for xx in range(256)])
    tau = np.array([256 * ((v // 256 + 1) % 3) + R[v % 256] for v in range(768)])
    return Adj, tau, R


def main():
    data = transport()
    (HERE.parent / 'data' / 'jordan_normal.json').write_text(json.dumps(data, indent=1) + '\n')
    Adj, tau, R = build_normal(data['z'], data['A'], data['B'])
    check_srg(Adj, (768, 312, 136, 120))
    print('normal-form graph: srg(768,312,136,120) verified on all pairs')
    print('tau = (i,x) -> (i+1, Rx) is an automorphism:', bool((Adj[np.ix_(tau, tau)] == Adj).all()),
          '; R fixed points:', int((R == np.arange(256)).sum()))
    from invariants768 import canonical_hash, to_igraph
    h1 = canonical_hash(to_igraph(Adj))
    h2 = canonical_hash(to_igraph(read_graph6(HERE.parent / 'graphs' / 'srg768_312_136_120_jordan_cayley.g6')))
    print('canonical hash normal form', h1, '; certificate graph', h2, '; isomorphic:', h1 == h2)
    print('Table 3 (sign of z at the vector with parameter u; w2 = w^2)')
    print('u: ' + ' '.join(f'({NAMES[p]},{NAMES[q]})' for p, q in U))
    z = np.array(data['z'])
    for k, S in SPACES.items():
        print(k, ' '.join('+' if z[x] > 0 else '-' for x in S))
    print('negative signs:', int((z < 0).sum()))
    for name in 'AB':
        print(name, '=', [[NAMES[e] for e in row] for row in data[name]])


if __name__ == '__main__':
    main()
