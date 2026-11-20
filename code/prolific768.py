"""Random members of the construction and their isomorphism classes.

python prolific768.py N [SEED]  -> JSON lines on stdout

A random member: three pairwise disjoint good four-sets of slopes (four points of a subline PG(1,4)),
a random one of the 12 C_5-invariant seeds on each, and a random translation c of the cross block B_0
(H of B_0 + c is 16 z_0 chi_c + 96 delta, so z_0 chi_c is again a seed, not C_5-invariant in general).
Each graph is verified on all pairs, then |Aut|, orbit sizes and a canonical hash are recorded.
"""
import itertools
import json
import random
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import srg768                                                  # noqa: E402
from verify_srg import check_srg                               # noqa: E402
from invariants768 import to_igraph, canonical_hash, automorphisms   # noqa: E402
from seed_search import SLOPES, coset_sums                     # noqa: E402
from srg768 import line, coset, N                              # noqa: E402


def all_seeds():
    signs = np.array(list(itertools.product((1, -1), repeat=12)))
    seeds = {}
    for four in itertools.combinations(SLOPES, 4):
        W = coset_sums(four) @ signs.T
        ok = ((W == 10) | (W == -6)).all(axis=0)
        if ok.any():
            seeds[four] = signs[ok]
    return seeds


def seed_vector(four, s):
    z = np.zeros(N, dtype=np.int64)
    for k, m in enumerate(four):
        pts = line(m)
        for a in range(1, 16):
            z[pts[a - 1]] = s[3 * k + coset(a)]
    return z


def main():
    n = int(sys.argv[1])
    rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
    seeds = all_seeds()
    fours = sorted(seeds)
    assert len(fours) == 340 and all(len(v) == 12 for v in seeds.values())
    # self-test: canonical hash invariant under relabeling
    A = srg768.cayley_graph()
    p = np.random.default_rng(1).permutation(len(A))
    assert canonical_hash(to_igraph(A)) == canonical_hash(to_igraph(A[np.ix_(p, p)]))
    for t in range(n):
        t0 = time.time()
        while True:
            chosen = rng.sample(fours, 3)
            if all(not set(a) & set(b) for a, b in itertools.combinations(chosen, 2)):
                break
        idx = [rng.randrange(12) for _ in range(3)]
        c = rng.randrange(N)
        zs = [seed_vector(f, seeds[f][i]) for f, i in zip(chosen, idx)]
        A = srg768.assemble(zs, translations=(c, 0, 0))
        check_srg(A, (768, 312, 136, 120))
        G = to_igraph(A)
        order, orbits, _ = automorphisms(G)
        print(json.dumps(dict(trial=t, fours=chosen, signings=idx, translation=c, aut=order,
                              orbits=orbits, hash=canonical_hash(G), seconds=round(time.time() - t0, 1))),
              flush=True)


if __name__ == '__main__':
    main()
