"""Invariants of the 768-vertex graphs: |Aut| and orbits (bliss), four-cliques through each vertex,
canonical hash, and isomorphism between the rebuilt graphs and the original files.

python invariants768.py [file.g6 ...]      (default: the two constructions of srg768.py and the files in graphs/)
"""
import hashlib
import sys
from collections import Counter
from pathlib import Path

import igraph as ig
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify_srg import read_graph6           # noqa: E402
import srg768                                  # noqa: E402


def to_igraph(A):
    ii, jj = np.nonzero(np.triu(A, 1))
    return ig.Graph(n=len(A), edges=list(zip(ii.tolist(), jj.tolist())))


def canonical_hash(G):
    perm = G.canonical_permutation(sh='fl')
    H = G.permute_vertices(perm)
    edges = sorted(tuple(sorted(e)) for e in H.get_edgelist())
    return hashlib.sha256(repr(edges).encode()).hexdigest()[:16]


def automorphisms(G):
    gens = G.automorphism_group(sh='fl')
    order = int(G.count_automorphisms(sh='fl'))
    parent = list(range(G.vcount()))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for g in gens:
        for i, j in enumerate(g):
            a, b = find(i), find(j)
            if a != b:
                parent[a] = b
    orbits = Counter(find(i) for i in range(G.vcount()))
    return order, sorted(orbits.values()), gens


def four_cliques_at(A, v):
    """Triangles in the neighbourhood of v = four-cliques containing v."""
    nb = np.flatnonzero(A[v])
    S = A[np.ix_(nb, nb)].astype(np.int64)
    return int(np.trace(S @ S @ S)) // 6


def report(name, A):
    G = to_igraph(A)
    order, orbits, _ = automorphisms(G)
    fc = [four_cliques_at(A, 256 * i) for i in range(3)]
    total = sum(four_cliques_at(A, v) for v in range(len(A))) // 4
    h = canonical_hash(G)
    print(f'{name}: |Aut| = {order}, orbit sizes {orbits}, four-cliques through vertex (i,0) for i=0,1,2: {fc}, '
          f'total four-cliques {total}, canonical hash {h}', flush=True)
    return h


if __name__ == '__main__':
    hashes = {}
    if len(sys.argv) > 1:
        for f in sys.argv[1:]:
            report(f, read_graph6(f))
        sys.exit()
    hashes['rebuilt shear'] = report('rebuilt three-orbit graph (srg768.shear_graph)', srg768.shear_graph())
    hashes['rebuilt cayley'] = report('rebuilt Cayley graph (srg768.cayley_graph)', srg768.cayley_graph())
    zs, _ = srg768.norm_cayley_seeds()
    hashes['rebuilt norm cayley'] = report('rebuilt norm Cayley graph (srg768.norm_cayley_seeds)', srg768.assemble(zs))
    hashes['rebuilt norm shear'] = report('rebuilt norm three-orbit graph (srg768.norm_threeorbit)',
                                          srg768.norm_threeorbit())
    for f in ['srg768_312_136_120_threeorbit.g6', 'srg768_312_136_120_cayley.g6',
              'srg768_312_136_120_norm_threeorbit.g6', 'srg768_312_136_120_norm_cayley.g6']:
        p = HERE.parent / 'graphs' / f
        if p.exists():
            hashes['graphs/' + f] = report('graphs/' + f, read_graph6(p))
    print('isomorphism classes by canonical hash:')
    for h in sorted(set(hashes.values())):
        print('  ', h, [k for k, v in hashes.items() if v == h])
