"""Discover four spread lines supporting a signed degree-96 Walsh block."""
import itertools
import json
import numpy as np
from construct768 import mul, sign_table


def projective_orbits():
    # Translations, scalar multiplications and inversion generate PGL(2,16).
    # Adjoin Frobenius for the semilinear action.
    inverse = [0] + [next(b for b in range(1,16) if mul(a,b)==1) for a in range(1,16)]
    maps = [[x^1 for x in range(16)]+[16],
            [mul(2,x) for x in range(16)]+[16],
            [16]+inverse[1:]+[0],
            [mul(x,x) for x in range(16)]+[16]]
    remaining = set(itertools.combinations(range(17),4))
    while remaining:
        first = min(remaining)
        remaining.remove(first)
        stack = [first]
        count = 0
        while stack:
            current = stack.pop()
            count += 1
            for perm in maps:
                other = tuple(sorted(perm[x] for x in current))
                if other in remaining:
                    remaining.remove(other)
                    stack.append(other)
        yield first, count


def signings(slopes, characters):
    orbit_sums = []
    for m in slopes:
        unused = set(range(1,16))
        while unused:
            x = min(unused)
            orbit = []
            for _ in range(5):
                unused.remove(x)
                orbit.append(x if m==16 else 16*x+mul(m,x))
                x = mul(8,x)
            orbit_sums.append(characters[:,orbit].sum(axis=1))
    table = np.array(orbit_sums).T
    total = table.sum(axis=1)
    good = []
    # At zero the Walsh transform must be 10. Each sign is repeated five times.
    for minus in itertools.combinations(range(12),5):
        values = total - 2*table[:,minus].sum(axis=1)
        if set(values) == {-6,10}:
            good.append(minus)
    return good


if __name__ == '__main__':
    characters = np.array(sign_table(),dtype=np.int64)
    for slopes, count in projective_orbits():
        print(json.dumps(dict(slopes=slopes, four_sets=count,
                              negative_orbit_indices=signings(slopes,characters))))
