"""Cayley members on C_2^8 x| C_3 with a 2-dimensional fixed space (the remaining group of this shape).

Every Cayley member of the generalized construction (Section 9) lives on V x| <R> with R of order 3. By the
fixed-space dimension of R these are C_2^8 x C_3 (dim 8, abelian: no PDS), dim 6 (impossible: a 4-dim support
subspace S meets the fixed space, so S and MS meet), dim 4 (the group G of Section 6), dim 0 (G_J of Section 9),
and dim 2: C_2^2 x (C_2^6 x| C_3). This script looks for the last.

Model: V = F_2^8 with the dot product (as in jordan768.py). M = I_2 + W + W + W, W = [[0,1],[1,1]] on the
coordinate pairs (2,3), (4,5), (6,7); M has order 3 and fixed space <e0, e1>.
We need four 4-dim subspaces A1..A4 with A_i cap M^j A_k = 0 for all (i,j,k) != (i,0,i), of type (b) or (c)
(type (a) carries no seeds), then a generalized seed on them (CP-SAT), then the graph.

python dim2_search.py SEED MAXTRIALS      (search; writes data/dim2_seed.json and checks/logs/dim2.g)
python dim2_search.py rebuild            (graph file and invariants from data/dim2_seed.json)
"""
import json
import random
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import jordan768 as J                       # noqa: E402

W = {0: 0, 1: 2, 2: 3, 3: 1}                 # x -> Wx on a pair: e_a -> e_b, e_b -> e_a + e_b


def M0(x):
    r = x & 3
    for k in (2, 4, 6):
        r |= W[(x >> k) & 3] << k
    return r


MMAP = [M0(x) for x in range(256)]
MIMAGES = [MMAP[1 << j] for j in range(8)]
assert all(MMAP[MMAP[MMAP[x]]] == x for x in range(256))
assert sum(MMAP[x] == x for x in range(256)) == 4      # fixed space of dimension 2


def mask(S):
    m = 0
    for x in S:
        m |= 1 << x
    return m


def images(S):
    S1 = [MMAP[x] for x in S]
    S2 = [MMAP[x] for x in S1]
    return mask(S), mask(S1), mask(S2)


def random_subspace(rng):
    while True:
        vecs = [rng.randrange(1, 256) for _ in range(4)]
        S = J.span(vecs)
        if len(S) == 16:
            return S


def search(seed, maxtrials):
    rng = random.Random(seed)
    trials = 0
    while trials < maxtrials:
        chosen, masks = [], []
        ok = True
        for k in range(4):
            for _ in range(2_000_000):
                trials += 1
                S = random_subspace(rng)
                m = images(S)
                if (m[0] & m[1]) != 1:
                    continue
                if all((a & b) == 1 for a in m for b in masks):
                    chosen.append(S)
                    masks.extend(m)
                    break
            else:
                ok = False
                break
        if not ok:
            continue
        typ = J.operator_type(*chosen)
        if typ in ('semisimple3', 'jordan'):
            return chosen, typ, trials
    return None, None, trials


def find_seed(support):
    from ortools.sat.python import cp_model
    pts = sorted(set(x for S in support for x in S) - {0})
    m = cp_model.CpModel()
    x = [m.NewBoolVar('') for _ in pts]
    y = [m.NewBoolVar('') for _ in range(256)]
    for v in range(256):
        co = [int(J.CHAR[v, b]) for b in pts]
        m.Add(sum(2 * c * xi for c, xi in zip(co, x)) - 16 * y[v] == sum(co) - 6)
    s = cp_model.CpSolver()
    s.parameters.num_workers = 4
    st = s.Solve(m)
    assert st in (cp_model.OPTIMAL, cp_model.FEASIBLE)
    z = np.zeros(256, dtype=np.int64)
    z[pts] = [2 * s.Value(t) - 1 for t in x]
    return z


def main():
    seed, maxtrials = int(sys.argv[1]), int(sys.argv[2])
    t = time.time()
    support, typ, trials = search(seed, maxtrials)
    print('trials', trials, 'time', round(time.time() - t, 1), 'type', typ, flush=True)
    if support is None:
        return
    z = find_seed([list(S) for S in support])
    print('seed found; Hz values', sorted(set((J.CHAR @ z).tolist())), flush=True)
    A, zs, sups, Cs, Mm = J.build(z, MIMAGES)
    from verify_srg import check_srg
    check_srg(A, (768, 312, 136, 120))
    print('graph: srg(768,312,136,120) on all pairs', flush=True)
    Minv = np.argsort(np.array(Mm))
    R = np.array([sum(int(J.DOT[x, Minv[1 << j]]) << j for j in range(8)) for x in range(256)])
    tau = np.array([256 * ((v // 256 + 1) % 3) + R[v % 256] for v in range(768)])
    print('tau automorphism:', bool((A[np.ix_(tau, tau)] == A).all()), '; R fixed points:',
          int((R == np.arange(256)).sum()), flush=True)
    from invariants768 import automorphisms, canonical_hash, four_cliques_at, to_igraph
    G = to_igraph(A)
    order, orbits, _ = automorphisms(G)
    print(f'|Aut| = {order}, orbits {orbits}, four-cliques at (i,0) {[four_cliques_at(A, 256 * i) for i in range(3)]}, '
          f'hash {canonical_hash(G)}', flush=True)
    out = Path(__file__).resolve().parent.parent / 'data' / 'dim2_seed.json'
    out.write_text(json.dumps({'M_basis_images': MIMAGES, 'support': [sorted(S) for S in support], 'type': typ,
                               'negative_frequencies': np.flatnonzero(z < 0).tolist(),
                               'positive_frequencies': np.flatnonzero(z > 0).tolist()}, indent=1) + '\n')
    # GAP: group generated by left multiplications, as in jordan768.py
    Rp = [np.arange(256), R, R[R]]

    def prod(g, h):
        return 256 * ((g // 256 + h // 256) % 3) + ((g % 256) ^ int(Rp[g // 256][h % 256]))
    gens = [[prod(h, g) for g in range(768)] for h in (1, 2, 4, 8, 16, 32, 64, 128, 256)]
    D = set(np.flatnonzero(A[0]).tolist())
    from collections import Counter

    def inverse(g):
        i, x = divmod(g, 256)
        j = (-i) % 3
        return 256 * j + int(Rp[j][x])
    cnt = Counter(prod(d, inverse(e)) for d in D for e in D)
    print('PDS counts ok:', cnt[0] == 312 and all(cnt[g] == (136 if g in D else 120) for g in range(1, 768)),
          '; left multiplications are automorphisms:',
          all((A[np.ix_(np.array(p), np.array(p))] == A).all() for p in gens), flush=True)
    gfile = Path(__file__).resolve().parent.parent / 'checks' / 'logs' / 'dim2.g'
    gfile.write_text('G := Group(' + ',\n'.join('PermList([' + ','.join(str(v + 1) for v in p) + '])' for p in gens)
                     + ');;\nPrint("|G| = ", Size(G), "  IdGroup ", IdGroup(G), "  |Z(G)| = ", Size(Centre(G)), "\\n");\n'
                     + 'Print("StructureDescription ", StructureDescription(G), "\\n");\nQUIT;\n')
    print('wrote', gfile.relative_to(gfile.parents[2]))


if __name__ == '__main__' and sys.argv[1:2] != ['rebuild']:
    main()


def graph6(A):
    n = len(A)
    bits = [int(A[i, j]) for j in range(n) for i in range(j)]
    bits += [0] * (-len(bits) % 6)
    out = [63 + 63, 63 + ((n >> 12) & 63), 63 + ((n >> 6) & 63), 63 + (n & 63)]
    out += [63 + int(''.join(map(str, bits[k:k + 6])), 2) for k in range(0, len(bits), 6)]
    return bytes(out).decode() + '\n'


def rebuild():
    """Rebuild the graph from data/dim2_seed.json, write graphs/srg768_312_136_120_dim2_cayley.g6, print invariants."""
    root = Path(__file__).resolve().parent.parent
    d = json.loads((root / 'data' / 'dim2_seed.json').read_text())
    z = np.zeros(256, dtype=np.int64)
    z[d['negative_frequencies']] = -1
    z[d['positive_frequencies']] = 1
    A, zs, sups, Cs, Mm = J.build(z, d['M_basis_images'])
    from verify_srg import check_srg, read_graph6
    check_srg(A, (768, 312, 136, 120))
    out = root / 'graphs' / 'srg768_312_136_120_dim2_cayley.g6'
    out.write_text(graph6(A))
    assert (read_graph6(out) == A).all()
    print('rebuilt from data/dim2_seed.json: srg on all pairs; wrote', out.name)
    print('support type (all 24 orders):', {J.operator_type(*[sups[0][k] for k in p])
                                             for p in __import__('itertools').permutations(range(4))})
    for i in range(3):
        n, types = J.jordan_invariant(Cs[i])
        print(f'fibre {i}: {n} four-dim subspaces in C_i + 0; four-sets by type {dict(types)}')
    from invariants768 import canonical_hash, to_igraph
    hj = canonical_hash(to_igraph(read_graph6(root / 'graphs' / 'srg768_312_136_120_jordan_cayley.g6')))
    print('canonical hash', canonical_hash(to_igraph(A)), '; Jordan graph of Section 9', hj)


if __name__ == '__main__' and sys.argv[1:2] == ['rebuild']:
    rebuild()
