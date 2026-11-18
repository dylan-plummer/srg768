"""Verify that graph6 files are strongly regular with given parameters (numpy, A^2 identity on all pairs).

python verify_srg.py v,k,lam,mu file.g6 [file.g6 ...]
Checks: symmetric 0/1, zero diagonal, A J = k J, and A^2 = k I + lam A + mu (J - I - A) entrywise.
"""
import sys
import numpy as np


def read_graph6(path):
    return decode_graph6(open(path, 'rb').read().strip())


def read_graph6_all(path):
    return [decode_graph6(line.strip()) for line in open(path, 'rb') if line.strip()]


def decode_graph6(s):
    if s.startswith(b'>>graph6<<'):
        s = s[10:]
    data = np.frombuffer(s, dtype=np.uint8).astype(np.int64) - 63
    assert (data >= 0).all() and (data < 64).all()
    if data[0] < 63:
        n, data = int(data[0]), data[1:]
    elif data[1] < 63:
        n = int(data[1] << 12 | data[2] << 6 | data[3]); data = data[4:]
    else:
        n = 0
        for x in data[2:8]:
            n = n << 6 | int(x)
        data = data[8:]
    bits = ((data[:, None] >> np.arange(5, -1, -1)) & 1).ravel()
    m = n * (n - 1) // 2
    assert len(bits) >= m and not bits[m:].any()
    A = np.zeros((n, n), dtype=np.uint8)
    # graph6: upper triangle in column order, i.e. (0,1),(0,2),(1,2),(0,3),...
    jj = np.repeat(np.arange(n), np.arange(n))
    ii = np.concatenate([np.arange(j) for j in range(n)]) if n > 1 else np.zeros(0, int)
    A[ii, jj] = bits[:m]
    A |= A.T
    return A


def check_srg(A, params, block=4096):
    v, k, lam, mu = params
    n = len(A)
    assert n == v, f'{n} vertices, expected {v}'
    assert (A == A.T).all() and not A.diagonal().any()
    assert (A.sum(1) == k).all(), 'not regular'
    Af = A.astype(np.float32)
    for start in range(0, n, block):
        S = Af[start:start + block] @ Af               # rows of A^2, exact in float32 (entries <= k < 2^24)
        Ab = A[start:start + block].astype(np.float32)
        expect = mu + (lam - mu) * Ab
        expect[np.arange(len(S)), np.arange(start, start + len(S))] = k
        assert np.array_equal(S, expect), f'SRG identity fails in rows {start}..'
    return True


if __name__ == '__main__':
    params = tuple(int(x) for x in sys.argv[1].split(','))
    for path in sys.argv[2:]:
        for t, A in enumerate(read_graph6_all(path)):
            check_srg(A, params)
            print(f'{path} [graph {t}]: srg{params} verified on all {len(A) * (len(A) - 1) // 2} pairs', flush=True)
