"""Check a partial difference set file (graphs/pds768_norm_cayley.txt) from its own header only.

python check_pds.py ../graphs/pds768_norm_cayley.txt [--gap out.g]

Parses R from the header, builds G = V x| <R> with (i,x)(j,y) = (i+j, x + R^i y), checks that D is inverse-closed,
avoids the identity and that every non-identity g is d e^{-1} exactly 136 times (g in D) or 120 times (g not in D).
Optionally writes a GAP script that identifies G.
"""
import re
import sys
from collections import Counter


def mul16(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a & 16:
            a ^= 0b10011
    return r


def main():
    path = sys.argv[1]
    text = open(path).read()
    m = re.search(r'R = \(\((\d+),(\d+)\),\((\d+),(\d+)\)\)', text)
    p, q, r, s = map(int, m.groups())

    def Rv(v):
        a, b = divmod(v, 16)
        return 16 * (mul16(p, a) ^ mul16(q, b)) + (mul16(r, a) ^ mul16(s, b))

    Rp = [list(range(256))]
    for _ in range(3):
        Rp.append([Rv(x) for x in Rp[-1]])
    assert Rp[3] == Rp[0] and Rp[1] != Rp[0], 'R must have order 3'

    def prod(g, h):
        i, x = divmod(g, 256)
        j, y = divmod(h, 256)
        return 256 * ((i + j) % 3) + (x ^ Rp[i][y])

    D = [256 * int(a) + int(b) for a, b in (ln.split() for ln in text.splitlines() if ln and not ln.startswith('#'))]
    Dset = set(D)
    assert len(D) == len(Dset) == 312 and 0 not in Dset
    inv = {}
    for g in range(768):
        inv[g] = next(h for h in range(768) if prod(g, h) == 0)
    assert all(inv[d] in Dset for d in D), 'D is not inverse-closed'
    cnt = Counter(prod(d, inv[e]) for d in D for e in D if d != e)
    assert all(cnt[g] == (136 if g in Dset else 120) for g in range(1, 768)), 'not a PDS'
    print(f'{path}: (768,312,136,120) partial difference set in V x| C_3 verified; '
          f'G abelian: {all(prod(g, h) == prod(h, g) for g in (1, 256) for h in (2, 512))}')
    if '--gap' in sys.argv:
        out = sys.argv[sys.argv.index('--gap') + 1]
        gens = [1, 2, 4, 8, 16, 32, 64, 128, 256]
        with open(out, 'w') as f:
            f.write('G := Group(' + ','.join('PermList([' + ','.join(str(prod(h, g) + 1) for g in range(768)) + '])'
                                            for h in gens) + ');;\n')
            f.write('Print("pds768_norm_cayley group: |G| = ", Size(G), " IdGroup ", IdGroup(G), " ", '
                    'StructureDescription(G), "\\n");\nQUIT;\n')


if __name__ == '__main__':
    main()
