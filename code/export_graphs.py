"""Regenerate the graph files in graphs/ that come from code/srg768.py (the norm-coordinate witnesses, a sample of
random family members, and the partial difference set of the norm Cayley graph).

The original witnesses (srg768_312_136_120_threeorbit.g6, srg768_312_136_120_cayley.g6 and the complement) come from
original/construct768.py and are not rewritten here.
"""
import random
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import srg768                                                # noqa: E402
from invariants768 import to_igraph, canonical_hash          # noqa: E402
from prolific768 import all_seeds, seed_vector               # noqa: E402

OUT = HERE.parent / 'graphs'


def graph6(A):
    n = len(A)
    assert n < 258048
    head = bytes([126, 63 + (n >> 12 & 63), 63 + (n >> 6 & 63), 63 + (n & 63)]) if n > 62 else bytes([63 + n])
    bits = [int(A[i, j]) for j in range(n) for i in range(j)]
    bits += [0] * (-len(bits) % 6)
    body = bytes(63 + int(''.join(map(str, bits[k:k + 6])), 2) for k in range(0, len(bits), 6))
    return head + body


def main():
    (OUT / 'srg768_312_136_120_norm_threeorbit.g6').write_bytes(graph6(srg768.norm_threeorbit()) + b'\n')
    zs, M = srg768.norm_cayley_seeds()
    A = srg768.assemble(zs)
    (OUT / 'srg768_312_136_120_norm_cayley.g6').write_bytes(graph6(A) + b'\n')

    # connection set of the norm Cayley graph: G = V x| <R>, R = M^{-T}, (i,x)(j,y) = (i+j, x + R^i y)
    (p, q), (r, s) = M
    det = srg768.mul(p, s) ^ srg768.mul(q, r)
    di = srg768.inv(det)
    R = ((srg768.mul(di, s), srg768.mul(di, r)), (srg768.mul(di, q), srg768.mul(di, p)))
    D = [int(g) for g in np.flatnonzero(A[0])]
    with open(OUT / 'pds768_norm_cayley.txt', 'w') as f:
        f.write('# Partial difference set (768,312,136,120) in G = V x| <R>, V = F_16^2 (vector (x,y) encoded 16x+y,\n'
                '# F_16 = F_2[t]/(t^4+t+1), integer bits = coefficients of 1,t,t^2,t^3), group law\n'
                f'# (i,x)(j,y) = (i+j mod 3, x + R^i y), R = (({R[0][0]},{R[0][1]}),({R[1][0]},{R[1][1]})) acting on columns.\n'
                '# One element per line: i x   (vertex 256 i + x of srg768_312_136_120_norm_cayley.g6)\n')
        for g in D:
            f.write(f'{g // 256} {g % 256}\n')

    seeds = all_seeds()
    fours = sorted(seeds)
    rng = random.Random(20261118)
    lines, hashes = [], set()
    while len(lines) < 20:
        chosen = rng.sample(fours, 3)
        if any(set(a) & set(b) for a in chosen for b in chosen if a != b):
            continue
        zs = [seed_vector(f, seeds[f][rng.randrange(12)]) for f in chosen]
        A = srg768.assemble(zs, translations=(rng.randrange(256), 0, 0))
        h = canonical_hash(to_igraph(A))
        if h in hashes:
            continue
        hashes.add(h)
        lines.append(graph6(A))
    (OUT / 'srg768_312_136_120_family_sample20.g6').write_bytes(b'\n'.join(lines) + b'\n')
    print('wrote', sorted(p.name for p in OUT.iterdir()))


if __name__ == '__main__':
    main()
