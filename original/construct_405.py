"""Rebuild and verify SRG(405,96,18,24) using only the standard library."""

import argparse
from collections import Counter
from itertools import combinations, product
from pathlib import Path


VERTICES = tuple(product(range(3), repeat=4))
WITHIN = (('0120', '1012'), ('0121', '1002'), ('0122', '1020'))
PORTS = (('0001', '0010', '0100'), ('0110', '1000', '1001'),
         ('0111', '1010', '1011'), ('0112', '1021', '1022'))


def values(u, basis):
    return tuple(sum(x * int(y) for x, y in zip(u, form)) % 3 for form in basis)


def construct():
    rows = [0] * 405
    profiles = [[values(u, port) for u in VERTICES] for port in PORTS]
    within = [[values(u, line) for line in WITHIN] for u in VERTICES]

    def add(x, y):
        rows[x] |= 1 << y
        rows[y] |= 1 << x

    for i in range(5):
        for u, v in combinations(range(81), 2):
            if any(a == b for a, b in zip(within[u], within[v])):
                add(i * 81 + u, i * 81 + v)
        for j in range(i + 1, 5):
            for u, (s, a, b) in enumerate(profiles[j - 1]):
                for v, (t, c, d) in enumerate(profiles[i]):
                    if (t == s and c == a) or ((t - s) % 3 == 1 and d == b):
                        add(i * 81 + u, j * 81 + v)
    return rows


def check(rows, parameters):
    v, k, lam, mu = parameters
    assert len(rows) == v and all(row.bit_count() == k for row in rows)
    assert all(not (row >> i & 1) and row >> v == 0 for i, row in enumerate(rows))
    counts = Counter()
    for j in range(v):
        for i in range(j):
            adjacent = rows[i] >> j & 1
            assert adjacent == rows[j] >> i & 1
            count = (rows[i] & rows[j]).bit_count()
            assert count == (lam if adjacent else mu), (i, j, count)
            counts[adjacent] += 1
    print(f'Verified SRG{parameters}. {counts[1]} adjacent pairs, {counts[0]} nonadjacent pairs.')


def graph6(rows):
    n = len(rows)
    if n < 63:
        prefix = chr(63 + n)
    elif n < 258048:
        prefix = '~' + ''.join(chr(63 + (n >> shift & 63)) for shift in (12, 6, 0))
    else:
        prefix = '~~' + ''.join(chr(63 + (n >> shift & 63)) for shift in (30, 24, 18, 12, 6, 0))
    bits = ''.join(str(rows[i] >> j & 1) for j in range(n) for i in range(j))
    bits += '0' * (-len(bits) % 6)
    return prefix + ''.join(chr(63 + int(bits[i:i + 6], 2)) for i in range(0, len(bits), 6)) + '\n'


def count_four_cliques(rows):
    count = 0
    for i, row in enumerate(rows):
        later = (row >> (i + 1)) << (i + 1)
        while later:
            next_vertex = later & -later
            later ^= next_vertex
            common = later & rows[next_vertex.bit_length() - 1]
            while common:
                third_vertex = common & -common
                common ^= third_vertex
                count += (common & rows[third_vertex.bit_length() - 1]).bit_count()
    return count


def check_clique_obstruction(rows):
    groups = ((333, 163, 324, 56, 62, 325, 112, 142), (162, 171, 301, 32, 109),
              (252, 243, 220, 244, 382))
    union = [v for group in groups for v in group]
    assert rows[0] >> 81 & 1
    assert len(set(union)) == len(union) and sum(1 << v for v in union) == rows[0] & rows[81]
    assert all(not (rows[v] >> w & 1) for group in groups for v, w in combinations(group, 2))
    print('Edge (0,81) lies in no clique of size 6 or greater.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path(__file__).parent)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    rows = construct()
    print(f'Four-cliques: {count_four_cliques(rows)}')
    check_clique_obstruction(rows)
    complement = [((1 << len(rows)) - 1) ^ row ^ (1 << i) for i, row in enumerate(rows)]
    for graph, parameters in [(rows, (405, 96, 18, 24)), (complement, (405, 308, 235, 231))]:
        check(graph, parameters)
        name = 'srg' + '_'.join(map(str, parameters)) + '.g6'
        (args.output / name).write_text(graph6(graph), encoding='ascii')


if __name__ == '__main__':
    main()
