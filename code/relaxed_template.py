"""Relaxed template on K_n: fibres an abelian group of order N (here elementary abelian p^m, m>=2), within-fibre
connection set C_i with Fourier values a0 (at 0), a1 on T_i (|T_i| = (n-1)s), a2 elsewhere; cross sets B_e with
|B^| = c on S_e (|S_e| = s), 0 elsewhere (beta != 0); supports of edges at a fibre disjoint, T_i = union of them.
Then eigenvalues a1 +- c on supports and a2 elsewhere: c = (t1-t2)/2, a1 = (t1+t2)/2, a2 = t2.
Solve trace/Parseval for (s, a0, kappa) for every open row of the local table."""
import csv, math
rows = []
for row in csv.DictReader(open(__import__('sys').argv[1] if len(__import__('sys').argv) > 1 else 'current_table.csv')):
    try: v, k, l, mu = int(row['v']), int(row['k']), int(row['lam']), int(row['mu'])
    except ValueError: continue
    rows.append((row['status'], v, k, l, mu, row['comments']))
def ppow(N):
    for p in range(2, N+1):
        if N % p == 0:
            e = 0
            while N % p == 0: N //= p; e += 1
            return (p, e) if N == 1 else None
for st, v, k, l, mu, com in rows:
    if st not in ('?',): continue
    D = (l-mu)**2 + 4*(k-mu); sq = math.isqrt(D)
    if sq*sq != D: continue
    t1, t2 = (l-mu+sq)//2, (l-mu-sq)//2
    if (t1 - t2) % 2: continue
    c = (t1-t2)//2; a1 = (t1+t2)//2; a2 = t2
    for n in range(2, v):
        if v % n: continue
        N = v // n; pe = ppow(N)
        if not pe or pe[1] < 2: continue
        for s in range(1, (N-1)//(n-1) + 1):
            # (1) trace: a0 = -( (n-1)s a1 + (N-1-(n-1)s) a2 )
            a0 = -((n-1)*s*a1 + (N-1-(n-1)*s)*a2)
            if not 0 <= a0 <= N-1: continue
            if a0*a0 + (n-1)*s*a1*a1 + (N-1-(n-1)*s)*a2*a2 != N*a0: continue
            if (k - a0) % (n-1) if n > 1 else True: continue
            kap = (k - a0)//(n-1)
            if not 0 < kap < N: continue
            if kap*kap + s*c*c != N*kap: continue
            if a0 - kap not in (t1, t2): continue
            print(f'srg({v},{k},{l},{mu}) eig {t1},{t2}: n={n} fibres of order {N}={pe[0]}^{pe[1]}, s={s}, '
                  f'within a0={a0} (two-char a1={a1} on {(n-1)*s}, a2={a2}), cross kappa={kap}, |B^|={c} on {s}  [{com[:40]}]')
