"""Check of Proposition 8.2 on line sums: for a prime p, the integer sequences n_0..n_{p-1} with 0 <= n_s <= p whose
Fourier transform has absolute value p at every t != 0 (equivalently, constant periodic autocorrelation
(kappa^2 - p^2)/p) are exactly the indicator of one point times p (a line coset) and its complement.
Exhaustive depth-first enumeration for p = 3, 5, 7.  python prime_obstruction.py"""

def flat_sequences(p):
    out = {}
    for kappa in range(p * p + 1):
        if (kappa * kappa - p * p) % p or (kappa * kappa + (p - 1) * p * p) % p:
            continue
        auto = (kappa * kappa - p * p) // p
        sq = (kappa * kappa + (p - 1) * p * p) // p
        n, sols = [0] * p, []

        def dfs(i, s, q):
            if q > sq or s > kappa:
                return
            if i == p:
                if s == kappa and q == sq and all(
                        sum(n[k] * n[(k + d) % p] for k in range(p)) == auto for d in range(1, p)):
                    sols.append(tuple(n))
                return
            for x in range(p + 1):
                if s + x + (p - i - 1) * p < kappa:
                    continue
                n[i] = x
                dfs(i + 1, s + x, q + x * x)
            n[i] = 0
        dfs(0, 0, 0)
        if sols:
            out[kappa] = sols
    return out


if __name__ == '__main__':
    for p in (3, 5, 7):
        res = flat_sequences(p)
        assert set(res) == {p, p * p - p}
        assert all(sorted(s) == [0] * (p - 1) + [p] for s in res[p])
        assert all(sorted(s) == [0] + [p] * (p - 1) for s in res[p * p - p])
        print(f'p={p}: flat line-sum sequences only for kappa = {p} (one full line) and {p*p-p} (its complement)')
