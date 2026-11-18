"""The Jordan-type Cayley graph srg(768,312,136,120) (data/jordan_seed.json), checked independently.

Model: V = F_2^8 (integers 0..255, bit j = coordinate j), pairing <x, y> = popcount(x & y) mod 2,
Walsh transform (Hf)(v) = sum_b f(b) (-1)^<v,b>.

python jordan768.py [--gap out.g]

1. The seed z (60 signed frequencies) is supported on four 4-dimensional subspaces of V meeting pairwise in 0,
   and Hz takes only the values 10 and -6. In coordinates where the first two subspaces are A + B = V and the
   third is the graph of phi: A -> B, the fourth is the graph of phi L for an operator L on A; L has minimal
   polynomial (X^2+X+1)^2 (not semisimple), so the four subspaces are not four lines of a Desarguesian spread.
2. M (8 basis images) has order 3, fixes only 0, and z, z o M^-1, z o M^-2 have disjoint supports. The graph of
   Construction 4.1 with these three seeds (and the general proof: the twelve support subspaces meet pairwise in 0)
   is built here, verified on all pairs, and compared with graphs/srg768_312_136_120_jordan_cayley.g6.
3. Cayley structure: tau (i,x) -> (i+1, R x), R = M^-T, is an automorphism; D = N(0,0) is a PDS in
   G = V x| <R>, group law (i,x)(j,y) = (i+j, x + R^i y). Optional GAP script for IdGroup and the centre.
4. |Aut| and orbits (bliss). Proposition 9.4 (not isomorphic to any graph of Construction 4.1): an isomorphism would
   be F_2-linear on fibres and map a fibre connection set C onto a set C' of Construction 4.1, so the number of
   4-dimensional subspaces in C + 0 and the types of their four-sets must agree. We compute them for each C_i of the
   Jordan graph (8 subspaces, 21 four-sets of type (c)); the paper's argument then needs nothing about C'. As a
   cross-check we also compute them for every C' of Construction 4.1 up to GammaL(2,16) (one subline four-set
   fixed, the other running over all disjoint ones): those with exactly 8 subspaces have no four-set of type (c).
"""
import json
import sys
from collections import Counter
from itertools import combinations, permutations
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify_srg import check_srg, read_graph6    # noqa: E402
import srg768                                     # noqa: E402

DOT = np.array([[bin(x & y).count('1') & 1 for y in range(256)] for x in range(256)], dtype=np.int64)
CHAR = 1 - 2 * DOT                                # CHAR[v, b] = (-1)^<v,b>


# ---------------------------------------------------------------- linear algebra over F_2 on 8-bit integers
def apply(images, x):
    r = 0
    for j, im in enumerate(images):
        if x >> j & 1:
            r ^= im
    return r


def span(vectors):
    S = {0}
    for v in vectors:
        S |= {s ^ v for s in S}
    return frozenset(S)


def annihilator(W):
    return frozenset(x for x in range(256) if all(not DOT[x, w] for w in W))


def subspaces_in(points, dim):
    """All subspaces of dimension dim contained in points | {0}."""
    P = set(points) | {0}
    level = {frozenset((0, p)) for p in P if p}
    for _ in range(dim - 1):
        nxt = set()
        for W in level:
            for p in P:
                if p not in W:
                    U = W | {w ^ p for w in W}
                    if U <= P:
                        nxt.add(frozenset(U))
        level = nxt
    return level


def basis_of(W):
    basis, S = [], {0}
    for w in sorted(W):
        if w not in S:
            basis.append(w)
            S |= {s ^ w for s in S}
    return basis


def solve_coords(basis):
    """Map x -> coordinate vector (int) of x in the given basis of V (8 vectors)."""
    table = {}
    for c in range(256):
        table[apply(basis, c)] = c
    assert len(table) == 256
    return table


def operator_type(A, B, C, D):
    """A, B, C, D: 4-dim subspaces meeting pairwise in 0. With V = A + B, C = graph(phi), D = graph(psi), return
    the type of L = phi^-1 psi on A: 'semisimple3' (L^2+L+1 = 0), 'quartic' (irreducible char. poly) or 'jordan'."""
    a, b = basis_of(A), basis_of(B)
    coords = solve_coords(a + b)

    def graph_map(Cs):
        m = {}
        for c in Cs:
            k = coords[c]
            m[k & 15] = k >> 4
        assert len(m) == 16 and len(set(m.values())) == 16
        return m
    phi, psi = graph_map(C), graph_map(D)
    phinv = {v: k for k, v in phi.items()}
    L = [phinv[psi[1 << j]] for j in range(4)]       # images of the basis of A, in A-coordinates

    def Lx(x):
        r = 0
        for j in range(4):
            if x >> j & 1:
                r ^= L[j]
        return r

    def q(x):                                         # (L^2 + L + 1) x
        return Lx(Lx(x)) ^ Lx(x) ^ x
    if all(q(x) == 0 for x in range(16)):
        return 'semisimple3'
    if all(q(q(x)) == 0 for x in range(16)):
        return 'jordan'
    # L and L + 1 are invertible, so the characteristic polynomial has no root in F_2; the only degree-4 options are
    # (X^2+X+1)^2 (handled above) and an irreducible quartic. Confirm irreducibility: no nonzero L-invariant
    # subspace of dimension 1 or 2, i.e. no x != 0 with x, Lx, L^2 x spanning less than dimension 3.
    for x in range(1, 16):
        y, w = Lx(x), Lx(Lx(x))
        assert len({0, x, y, x ^ y, w, w ^ x, w ^ y, w ^ x ^ y}) == 8, 'reducible but not of type (X^2+X+1)^2'
    return 'quartic'


def jordan_invariant(C):
    """(number of 4-subspaces inside C + 0, Counter of operator types over ordered 4-tuples meeting pairwise in 0)."""
    subs = sorted(subspaces_in(C, 4), key=sorted)
    types = Counter()
    for quad in combinations(subs, 4):
        if all(len(X & Y) == 1 for X, Y in combinations(quad, 2)):
            ts = {operator_type(*[quad[k] for k in p]) for p in permutations(range(4))}
            assert len(ts) == 1, 'operator type depends on the order'
            types[ts.pop()] += 1
    return len(subs), types


# ---------------------------------------------------------------- the graph
def load():
    d = json.loads((HERE.parent / 'data' / 'jordan_seed.json').read_text())
    z = np.zeros(256, dtype=np.int64)
    z[d['negative_frequencies']] = -1
    z[d['positive_frequencies']] = 1
    return z, d['frequency_basis_images']


def support_subspaces(z):
    S = set(np.flatnonzero(z).tolist())
    subs = [W for W in subspaces_in(S, 4)]
    return subs


def build(z, M):
    Mmap = [apply(M, x) for x in range(256)]
    zs = [z]
    for _ in range(2):
        nz = np.zeros(256, dtype=np.int64)
        for b in np.flatnonzero(zs[-1]):
            nz[Mmap[b]] = zs[-1][b]                   # z_{i+1} = z_i o M^-1
        zs.append(nz)
    sup = [support_subspaces(zi) for zi in zs]
    allsubs = [W for s in sup for W in s]
    assert all(len(s) == 4 for s in sup)
    assert all(len(X & Y) == 1 for X, Y in combinations(allsubs, 2)), 'twelve subspaces must meet pairwise in 0'
    Bs = []
    for zi in zs:
        Hz = CHAR @ zi
        assert set(Hz.tolist()) == {10, -6}
        Bs.append(np.flatnonzero(Hz == 10))
    Cs = []
    for i in range(3):
        C = set()
        for W in sup[i - 1] + sup[i]:
            C |= annihilator(W) - {0}
        assert len(C) == 120
        Cs.append(sorted(C))
    A = np.zeros((768, 768), dtype=np.uint8)
    x = np.arange(256)
    for i in range(3):
        for c in Cs[i]:
            A[256 * i + x, 256 * i + (x ^ c)] = 1
        j = (i + 1) % 3
        for b in Bs[i]:
            A[256 * i + x, 256 * j + (x ^ b)] = 1
            A[256 * j + (x ^ b), 256 * i + x] = 1
    return A, zs, sup, Cs, Mmap


def main():
    z, M = load()
    sup = support_subspaces(z)
    print('support: 4-dim subspaces in supp z + 0:', len(sup), '; meet pairwise in 0:',
          all(len(X & Y) == 1 for X, Y in combinations(sup, 2)), '; |supp z| =', int(np.count_nonzero(z)))
    print('Hz values:', dict(Counter((CHAR @ z).tolist())))
    print('operator type of the support (in the stored order of subspaces, all 24 orders):',
          Counter(operator_type(*[sup[k] for k in p]) for p in permutations(range(4))))
    A, zs, sups, Cs, Mmap = build(z, M)
    Mmap = np.array(Mmap)
    print('M: order 3', bool((Mmap[Mmap[Mmap]] == np.arange(256)).all()), '; fixed points', int((Mmap == np.arange(256)).sum()))
    check_srg(A, (768, 312, 136, 120))
    check_srg(1 - A - np.eye(768, dtype=np.uint8), (768, 455, 262, 280))
    print('rebuilt graph: srg(768,312,136,120) on all pairs; complement srg(768,455,262,280)')
    F = read_graph6(HERE.parent / 'graphs' / 'srg768_312_136_120_jordan_cayley.g6')
    print('identical to graphs/srg768_312_136_120_jordan_cayley.g6 (same labelling):', bool((F == A).all()))

    # Cayley structure. R = M^-T: <R x, y> = <x, M^-1 y>.
    Minv = np.argsort(Mmap)
    R = np.array([sum((int(DOT[x, Minv[1 << j]]) << j) for j in range(8)) for x in range(256)])
    assert all(DOT[R[x], y] == DOT[x, Minv[y]] for x in range(256) for y in range(0, 256, 7))
    Rp = [np.arange(256), R, R[R]]
    tau = np.array([256 * ((v // 256 + 1) % 3) + R[v % 256] for v in range(768)])
    print('tau is an automorphism:', bool((A[np.ix_(tau, tau)] == A).all()))

    def prod(g, h):
        return 256 * ((g // 256 + h // 256) % 3) + ((g % 256) ^ int(Rp[g // 256][h % 256]))

    def inverse(g):
        i, x = divmod(g, 256)
        j = (-i) % 3
        return 256 * j + int(Rp[j][x])
    D = np.flatnonzero(A[0]).tolist()
    Dset = set(D)
    assert all(inverse(d) in Dset for d in D) and 0 not in Dset
    cnt = Counter(prod(d, inverse(e)) for d in D for e in D)
    ok = cnt[0] == 312 and all(cnt[g] == (136 if g in Dset else 120) for g in range(1, 768))
    # left translations by G are automorphisms: g -> h g
    lt = all((A[np.ix_(p, p)] == A).all() for p in
             [np.array([prod(h, g) for g in range(768)]) for h in (1, 2, 4, 8, 16, 32, 64, 128, 256)])
    print('D = N(0): inverse-closed, PDS counts (312 at 1, 136 on D, 120 off D):', ok, '; left multiplication by generators of G are automorphisms and Cay(G,D) = graph:', lt)

    if '--gap' in sys.argv:
        out = Path(sys.argv[sys.argv.index('--gap') + 1])

        def gap_perm(p):
            return 'PermList([' + ','.join(str(int(v) + 1) for v in p) + '])'
        gens = [np.array([prod(h, g) for g in range(768)]) for h in (1, 2, 4, 8, 16, 32, 64, 128, 256)]
        Dp = [d + 1 for d in D]
        lines = ['G := Group(' + ',\n'.join(gap_perm(p) for p in gens) + ');;',
                 'Print("|G| = ", Size(G), "  IdGroup ", IdGroup(G), "  |Z(G)| = ", Size(Centre(G)), '
                 '"  IsAbelian ", IsAbelian(G), "\\n");',
                 'Print("StructureDescription ", StructureDescription(G), "\\n");',
                 'QUIT;']
        out.write_text('\n'.join(lines) + '\n')
        print('wrote', out)

    # Automorphisms.
    from invariants768 import automorphisms, canonical_hash, four_cliques_at, to_igraph
    G = to_igraph(A)
    order, orbits, _ = automorphisms(G)
    print(f'|Aut| = {order}, orbits {orbits}, four-cliques through (i,0): '
          f'{[four_cliques_at(A, 256 * i) for i in range(3)]}, canonical hash {canonical_hash(G)}', flush=True)

    # GL(8,2)-invariant of fibre connection sets.
    for i in range(3):
        n, types = jordan_invariant(Cs[i])
        print(f'Jordan graph, fibre {i}: {n} four-dim subspaces in C_i + 0; four-sets meeting pairwise in 0 by operator type (the same for all 24 orders) {dict(types)}', flush=True)
    base = (0, 1, 6, 7)                               # the subline {0, 1, w, w^2}
    seen = Counter()
    for four in combinations(range(17), 4):
        if set(four) & set(base):
            continue
        if srg768_cross_class(four) != (6, 7):
            continue
        C = set()
        for m in base + four:
            C |= set(srg768.line(srg768.annihilator_slope(m)))
        n, types = jordan_invariant(C)
        seen[(n, types.get('jordan', 0) > 0)] += 1
    print('Construction 4.1 fibre sets C\' (T = {0,1,w,w^2} u S, S running over the', sum(seen.values()),
          'subline four-sets disjoint from it): (number of 4-dim subspaces in C\' + 0, some Jordan 4-tuple) -> count')
    for k, v in sorted(seen.items()):
        print(f'   {k}: {v}')
    print('a set C\' with exactly 8 subspaces and a Jordan 4-tuple exists:', any(k == (8, True) for k in seen))


def srg768_cross_class(four):
    from seed_search import cross_ratio_class
    return cross_ratio_class(list(four))


if __name__ == '__main__':
    main()
