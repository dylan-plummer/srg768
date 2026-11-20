"""Construct and exactly verify the 768-vertex spread and Cayley witnesses."""
import argparse
from collections import Counter
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'construction_frontier_20261117'))
from construct_405 import check, graph6

PARAMETERS = (768, 312, 136, 120)
COMPLEMENT = (768, 455, 262, 280)
SLOPES = (0, 1, 2, 10)
NEGATIVE = {0, 1, 3, 6, 11}


def mul(a, b):
    """Polynomial-basis F16 using X^4 + X + 1."""
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 16:
            a ^= 19
    return result


def trace(a):
    value = 0
    for _ in range(4):
        value ^= a
        a = mul(a, a)
    assert value in (0, 1)
    return value


def sign_table():
    tr = [[trace(mul(i, j)) for j in range(16)] for i in range(16)]
    return [[1-2*(tr[x>>4][y>>4] ^ tr[x&15][y&15]) for y in range(256)] for x in range(256)]


def seed():
    result = {}
    for m, slope in enumerate(SLOPES):
        remaining = set(range(1, 16))
        for j in range(3):
            x = min(remaining)
            for _ in range(5):
                remaining.remove(x)
                result[16*x + mul(slope, x)] = -1 if 3*m+j in NEGATIVE else 1
                x = mul(8, x)
        assert not remaining
    return result


def frequency_step(v):
    x, y = divmod(v, 16)
    return 16*x + (mul(5, x) ^ mul(6, y))


def vertex_step(v):
    # Inverse transpose of [[1,0],[5,6]] for the field trace pairing.
    x, y = divmod(v, 16)
    z = mul(7, y)
    return 16*(x ^ mul(5, z)) + z


def spectral_edges(kind):
    first = seed()
    if kind == 'cayley':
        second = {frequency_step(v): s for v, s in first.items()}
        third = {frequency_step(v): s for v, s in second.items()}
        assert {frequency_step(v): s for v, s in third.items()} == first
        return [first, second, third]
    def shear(v, h):
        x, y = divmod(v, 16)
        return 16*x + (y ^ mul(h, x))
    return [{shear(v, h): s for v, s in first.items()} for h in (0, 4, 13)]


def construct(kind):
    chars = sign_table()
    edges = spectral_edges(kind)
    assert len(set.union(*(set(e) for e in edges))) == 180
    within, across = [], []
    for i, edge in enumerate(edges):
        transform = [sum(row[v]*s for v, s in edge.items()) for row in chars]
        assert Counter(transform) == {-6:160, 10:96}
        across.append({j for j, w in enumerate(transform) if w == 10})
        support = set(edge) | set(edges[i-1])
        values = [sum(row[v] for v in support) for row in chars]
        assert Counter(values) == {120:1, 8:120, -8:135}
        within.append({j for j, w in enumerate(values) if w == 8})
    rows = []
    for i in range(3):
        for x in range(256):
            neighbours = [(i, within[i]), ((i+1)%3, across[i]), ((i-1)%3, across[i-1])]
            rows.append(sum(1 << (256*j + (x ^ d)) for j, differences in neighbours for d in differences))
    return rows


def product(g, h):
    i, x = divmod(g, 256)
    j, y = divmod(h, 256)
    for _ in range(i):
        y = vertex_step(y)
    return 256*((i+j)%3) + (x ^ y)


def verify_cayley(rows):
    assert all(vertex_step(vertex_step(vertex_step(x))) == x for x in range(256))
    assert all(vertex_step(x ^ (1<<j)) == vertex_step(x) ^ vertex_step(1<<j)
               for x in range(256) for j in range(8))
    neighbours = [g for g in range(768) if rows[0] >> g & 1]
    for g, row in enumerate(rows):
        assert row == sum(1 << product(g, h) for h in neighbours)
    inverses = {g: next(h for h in neighbours if product(g, h)==product(h, g)==0)
                for g in neighbours}
    counts = Counter(product(g, inverses[h]) for g in neighbours for h in neighbours)
    assert counts[0] == 312
    assert all(counts[g] == (136 if rows[0]>>g & 1 else 120) for g in range(1, 768))
    return len(neighbours)


def local_four_cliques(rows, vertex):
    """Count triangles in the neighborhood, each in increasing vertex order."""
    candidates = rows[vertex]
    count = 0
    while candidates:
        first = candidates & -candidates
        candidates ^= first
        later = candidates & rows[first.bit_length()-1]
        while later:
            second = later & -later
            later ^= second
            count += (later & rows[second.bit_length()-1]).bit_count()
    return count


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, default=HERE)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    for kind in ['shear', 'cayley']:
        rows = construct(kind)
        check(rows, PARAMETERS)
        signature = [local_four_cliques(rows, 256*i) for i in range(3)]
        print(kind, 'four-cliques through a representative of each fiber', signature)
        if kind == 'cayley':
            print('Cayley connection set verified, size', verify_cayley(rows))
        name = 'srg768_312_136_120' + ('_cayley' if kind=='cayley' else '')
        (args.output/(name+'.g6')).write_text(graph6(rows), encoding='ascii')
        if kind == 'shear':
            comp = [((1 << 768)-1) ^ row ^ (1 << i) for i, row in enumerate(rows)]
            check(comp, COMPLEMENT)
            (args.output/'srg768_455_262_280.g6').write_text(graph6(comp), encoding='ascii')


if __name__ == '__main__':
    main()
