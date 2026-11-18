"""Parameter analysis: which (n, q, r) can the construction support?

Generalisation: n fibres, each a copy of F_q^2, joined along all pairs (fibre graph K_n). Edge e carries a seed
z_e: supported on r punctured spread lines with |z_e| = 1 there, and a set B_e with H 1_{B_e} = q z_e + k delta_0.
Edges at a common fibre have disjoint supports; fibre i uses the annihilators of the (n-1) r lines of its edges.

Block on a character beta in S_e:   [[q - (n-1) r, q z], [q conj z, q - (n-1) r]]  -> 2q - (n-1) r, -(n-1) r
Character in no support at fibre i: -(n-1) r.
Constant characters: (n-1) r (q-1) I + k (J - I): nonprincipal eigenvalue (n-1) r (q-1) - k.
Parseval for B_e: k^2 - q^2 k + r (q-1) q^2 = 0.

So we need integer k, (n-1) r (q-1) - k in {2q - (n-1) r, -(n-1) r}, and r * chi'(K_n) <= q + 1 lines.
Prints every solution with q a prime power <= QMAX and the resulting srg parameters, with the status in the local
copy of Brouwer's table if available.
"""
import csv
import math
import sys
from pathlib import Path

QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 512


def prime_power(q):
    for p in range(2, q + 1):
        if q % p == 0:
            while q % p == 0:
                q //= p
            return q == 1
    return False


def table():
    path = Path(__file__).resolve().parents[2] / 'current_table.csv'
    out = {}
    if path.exists():
        for row in csv.DictReader(open(path)):
            try:
                out[(int(row['v']), int(row['k']), int(row['lam']), int(row['mu']))] = (row['status'], row['comments'])
            except ValueError:
                pass
    return out


def main():
    tab = table()
    print('n  q   r  k_cross  srg parameters            eigenvalues   lines used/available   table')
    for q in range(2, QMAX + 1):
        if not prime_power(q):
            continue
        for n in range(2, 4 * q):
            d = n - 1
            colours = d if n % 2 == 0 else n
            if n == 2:
                colours = 1
            for r in range(1, (q + 1) // colours + 1):
                disc = q ** 4 - 4 * r * (q - 1) * q * q
                if disc < 0:
                    continue
                s = math.isqrt(disc)
                if s * s != disc:
                    continue
                for k in {(q * q + s) // 2, (q * q - s) // 2}:
                    if (q * q + s) % 2 or not (0 < k < q * q):
                        continue
                    th1, th2 = 2 * q - d * r, -d * r
                    if d * r * (q - 1) - k not in (th1, th2):
                        continue
                    v = n * q * q
                    K = d * r * (q - 1) + d * k
                    mu = K + th1 * th2
                    lam = mu + th1 + th2
                    if not (0 < K < v - 1):
                        continue
                    status = tab.get((v, K, lam, mu)) or tab.get((v, v - K - 1, v - 2 - 2 * K + mu, v - 2 * K + lam))
                    print(f'{n:<2} {q:<3} {r:<2} {k:<8} srg({v},{K},{lam},{mu})'.ljust(46)
                          + f'{th1},{th2}'.ljust(14) + f'{colours * r}/{q + 1}'.ljust(23)
                          + (f'{status[0]} {status[1][:60]}' if status else ('v>1300' if v > 1300 else '')))


if __name__ == '__main__':
    main()
