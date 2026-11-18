"""Exhaustive search for C_5-invariant seeds.

A seed is z: V -> {0, +-1} supported on four punctured spread lines with Walsh transform in {10, -6}.
Here z is required to be constant on the three C_5-cosets of each line (C_5 = <t^3> < F_16^* acting by scalars),
so a candidate is (four slopes, 12 signs). We test all 2380 four-sets and all 2^12 sign vectors.

Also records the action of the semilinear group PGammaL(2,16) on four-sets of slopes (two orbits, by cross-ratio)
and prints the twelve signings on the base slopes (0, 1, t, t^3+t) as a table sign[m][coset].
"""
import itertools
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from srg768 import CHI, INF, LOG, coset, inv, line, mul   # noqa: E402

SLOPES = list(range(16)) + [INF]


def cross_ratio_class(four):
    """Orbit invariant of PGammaL(2,16) on 4-subsets of PG(1,16): the set of cross-ratios is a {lambda, 1/lambda,
    1+lambda, ...} orbit; we compute the 6-set of cross-ratios over all orderings and take its Frobenius orbit."""
    def cr(a, b, c, d):
        # cross-ratio (a,b;c,d) = (c-a)(d-b)/((c-b)(d-a)) with infinity handled by dropping factors
        def diff(x, y):
            return None if INF in (x, y) else x ^ y
        num = [diff(c, a), diff(d, b)]
        den = [diff(c, b), diff(d, a)]
        val = 1
        for f in num:
            if f is not None:
                val = mul(val, f)
        for f in den:
            if f is not None:
                val = mul(val, inv(f))
        return val
    vals = frozenset(cr(*p) for p in itertools.permutations(four))
    frob = [vals]
    for _ in range(3):
        frob.append(frozenset(mul(x, x) for x in frob[-1]))
    return min(tuple(sorted(f)) for f in frob)


def coset_sums(four):
    """12 columns: for slope k and coset j, the vector beta -> sum over a in coset j of chi_{(a, m a)}(beta)."""
    cols = []
    for m in four:
        pts = line(m)
        for j in range(3):
            idx = [pts[a - 1] for a in range(1, 16) if coset(a) == j]
            cols.append(CHI[:, idx].sum(axis=1))
    return np.array(cols).T                       # 256 x 12


def main():
    signs = np.array(list(itertools.product((1, -1), repeat=12)))      # 4096 x 12
    by_class = Counter()
    good_by_class = Counter()
    sizes = Counter()
    base = (0, 1, 2, 10)
    base_good = None
    for four in itertools.combinations(SLOPES, 4):
        cls = cross_ratio_class(four)
        by_class[cls] += 1
        W = coset_sums(four) @ signs.T                                  # 256 x 4096
        ok = ((W == 10) | (W == -6)).all(axis=0)
        good_by_class[cls] += int(ok.sum())
        sizes[(cls, int(ok.sum()))] += 1
        if four == base:
            base_good = signs[ok]
    print('four-sets of slopes by cross-ratio class (PGammaL(2,16) orbit invariant):')
    for cls, n in by_class.items():
        print(f'  class {cls}: {n} four-sets, {good_by_class[cls]} seeds in total')
    print('distribution (class, seeds per four-set) -> number of four-sets:', dict(sizes))
    print(f'base slopes {base}: {len(base_good)} seeds; sign[m][coset j], coset j = log_t(a) mod 3:')
    for s in base_good:
        print('   ', ' | '.join(f'm={m}: ' + ' '.join('+' if x > 0 else '-' for x in s[3 * k:3 * k + 3])
                                 for k, m in enumerate(base)))
    with open(Path(__file__).resolve().parent.parent / 'checks' / 'logs' / 'seeds_base.json', 'w') as f:
        json.dump({'slopes': base, 'coset_index': 'log_t(a) mod 3',
                   'signings': base_good.tolist()}, f)


if __name__ == '__main__' and 'gl' not in sys.argv:
    main()


def gl_orbit_check():
    """Every C_5-invariant seed is the image of the norm seed z(a, mu a) = eps_mu (-1)^[N(a)=1] (mu in F_4,
    eps = -1 only at mu = 0) under some M in GL(2,16) acting on frequencies."""
    from srg768 import norm_seed, transport, F4
    z = norm_seed(0, 0, 1)
    pts = np.flatnonzero(z)
    images = set()
    for p, q, r, s in itertools.product(range(16), repeat=4):
        if mul(p, s) ^ mul(q, r) == 0:
            continue
        w = transport(z, ((p, q), (r, s)))
        images.add(tuple(w.tolist()))
    signs = np.array(list(itertools.product((1, -1), repeat=12)))
    found = set()
    for four in itertools.combinations(SLOPES, 4):
        W = coset_sums(four) @ signs.T
        for col in np.flatnonzero(((W == 10) | (W == -6)).all(axis=0)):
            zz = np.zeros(256, dtype=np.int64)
            for k, m in enumerate(four):
                for a in range(1, 16):
                    zz[line(m)[a - 1]] = signs[col][3 * k + coset(a)]
            found.add(tuple(zz.tolist()))
    print(f'GL(2,16)-images of the norm seed: {len(images)}; C_5-invariant seeds from the exhaustive search: '
          f'{len(found)}; equal: {images == found}')


if __name__ == '__main__' and 'gl' in sys.argv:
    gl_orbit_check()
