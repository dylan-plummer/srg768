r"""The Cayley graph as a partial difference set in G = V x|_R C_3, and its automorphism group in GAP.

1. Builds G with product (i, x)(j, y) = (i + j, x + R^i y), where R = M^{-T} and M(a, b) = (a, 5a + 6b) acts on
   frequencies. Checks that the spectral construction srg768.cayley_graph() equals Cay(G, D), D = N(identity),
   and that D is a regular (768, 312, 136, 120) partial difference set: every element of D \ {1} is written
   136 times as d e^{-1}, every other non-identity element 120 times.
2. Writes a GAP script (checks/logs/cayley768.g) that
   - identifies G (IdGroup) and its structure,
   - computes Aut(Gamma) from bliss generators, its order, solvability, and its regular subgroups of order 768.
"""
import sys
from collections import Counter
from pathlib import Path

import igraph as ig
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import srg768                                       # noqa: E402
from srg768 import apply_matrix, vertex_map, M_FREQ, N  # noqa: E402

R = vertex_map()


def R_pow(i, y):
    for _ in range(i % 3):
        y = apply_matrix(R, y)
    return y


def prod(g, h):
    i, x = divmod(g, N)
    j, y = divmod(h, N)
    return N * ((i + j) % 3) + (x ^ R_pow(i, y))


def main():
    # R is linear of order 3 and preserves the trace pairing with M
    assert all(R_pow(3, x) == x for x in range(N))
    for v in range(0, N, 7):
        for beta in range(0, N, 5):
            assert srg768.CHI[apply_matrix(M_FREQ, beta), apply_matrix(R, v)] == srg768.CHI[beta, v]
    A = srg768.cayley_graph()
    D = np.flatnonzero(A[0])
    assert len(D) == 312
    inverse = {}
    for g in range(3 * N):
        for h in range(3 * N):
            if prod(g, h) == 0:
                inverse[g] = h
                break
    table = np.array([[prod(g, h) for h in range(3 * N)] for g in range(3 * N)], dtype=np.int32)
    # Cay(G, D): g ~ g d
    C = np.zeros_like(A)
    for g in range(3 * N):
        C[g, table[g, D]] = 1
    assert (C == A).all(), 'spectral graph differs from the Cayley graph'
    assert all(inverse[d] in set(D.tolist()) for d in D)
    cnt = Counter(int(table[d, inverse[e]]) for d in D for e in D)
    inD = set(D.tolist())
    assert cnt[0] == 312
    assert all(cnt[g] == (136 if g in inD else 120) for g in range(1, 3 * N))
    print('Cay(G, D) equals the spectral construction; D is a (768,312,136,120) PDS in G = V x| C_3')
    print('G is nonabelian:', any(table[g, h] != table[h, g] for g in (1, 256, 257) for h in (256, 3, 513)))

    # GAP input: right regular action g -> g h of the generators, and Aut(Gamma) generators from bliss
    G = ig.Graph(n=3 * N, edges=[tuple(e) for e in np.argwhere(np.triu(A, 1)).tolist()])
    aut_gens = G.automorphism_group(sh='fl')
    gens = [1, 2, 4, 8, 16, 32, 64, 128, 256]                     # translations and the C_3 generator
    def gap_perm(p):
        return 'PermList([' + ','.join(str(x + 1) for x in p) + '])'
    left = [[int(table[h, g]) for g in range(3 * N)] for h in gens]   # left multiplication: automorphisms of Cay(G,D)
    out = HERE.parent / 'checks' / 'logs' / 'cayley768.g'
    with open(out, 'w') as f:
        f.write('LoadPackage("smallgrp");;\n')
        f.write('G := Group(' + ','.join(gap_perm(p) for p in left) + ');;\n')
        f.write('Print("|G| = ", Size(G), "  IsAbelian ", IsAbelian(G), "\\n");\n')
        f.write('Print("IdGroup(G) = ", IdGroup(G), "  StructureDescription ", StructureDescription(G), "\\n");\n')
        f.write('Aut := Group(' + ','.join(gap_perm(p) for p in aut_gens) + ');;\n')
        f.write('Print("|Aut| = ", Size(Aut), "  solvable ", IsSolvable(Aut), "  G <= Aut ", IsSubgroup(Aut, G), "\\n");\n')
        f.write('Print("Aut structure ", StructureDescription(Aut), "\\n");\n')
        f.write('Print("vertex stabiliser ", StructureDescription(Stabilizer(Aut, 1)), "\\n");\n')
        f.write('H := HallSubgroup(Aut, [2,3]);;\n')
        f.write('Print("Hall {2,3}-subgroup order ", Size(H), " regular ", IsRegular(H, [1..768]), '
                '" Id ", IdGroup(H), "\\n");\n')
        f.write('cls := Filtered(ConjugacyClassesSubgroups(Aut), c -> Size(Representative(c)) = 768);;\n')
        f.write('Print("classes of subgroups of order 768: ", Length(cls), " regular: ", '
                'List(cls, c -> IsRegular(Representative(c), [1..768])), " ids ", '
                'List(cls, c -> IdGroup(Representative(c))), "\\n");\n')
        f.write('QUIT;\n')
    print('wrote', out)


if __name__ == '__main__':
    main()
