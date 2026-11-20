"""Strongly regular graphs srg(768,312,136,120) from three copies of F_16^2.

Written from the mathematical description only (no code shared with original/).

Conventions
-----------
F = F_2[t]/(t^4+t+1). An integer 0..15 encodes the element sum_i bit_i t^i.
V = F^2. The vector (x, y) is encoded as 16*x + y.
Characters: chi_(a,b)(x,y) = (-1)^Tr(a x + b y), where Tr is the absolute trace F -> F_2.
Walsh transform: (H f)(beta) = sum_alpha f(alpha) chi_alpha(beta), so H^2 = 256 I.

Spread lines: L_m = {(a, m a)} for m in F, and L_inf = {(0, a)}. Line L_m is identified with its slope m (16 = inf).

A *seed* is a function z: V -> {0, +1, -1} supported on four punctured spread lines with
  H z in {10, -6} pointwise.
Then b = (H z + 6)/16 is the indicator of a set B of size 96 with H 1_B = 16 z + 96 delta_0.

Assembly (three fibers i in Z/3, vertex (i, x)):
  * edge {i, i+1} carries a seed z_i with support S_i; the supports are pairwise disjoint;
  * (i, x) ~ (i+1, y) iff x + y in B_i;
  * inside fiber i: T_i = S_{i-1} u S_i (eight punctured lines) and (i, x) ~ (i, y) iff x + y lies in
    one of the eight punctured lines L^perp, L in T_i.
"""
import itertools
import numpy as np

# ---------------------------------------------------------------- F_16
EXP = [0] * 30
LOG = [None] * 16
_x = 1
for _i in range(15):
    EXP[_i] = EXP[_i + 15] = _x
    LOG[_x] = _i
    _x <<= 1
    if _x & 16:
        _x ^= 0b10011          # t^4 = t + 1
assert sorted(EXP[:15]) == list(range(1, 16))


def mul(a, b):
    return 0 if a == 0 or b == 0 else EXP[LOG[a] + LOG[b]]


def inv(a):
    assert a
    return EXP[(15 - LOG[a]) % 15]


def tr(a):
    """Absolute trace a + a^2 + a^4 + a^8."""
    s, p = 0, a
    for _ in range(4):
        s ^= p
        p = mul(p, p)
    assert s in (0, 1)
    return s


N = 256
INF = 16
TRMUL = np.array([[tr(mul(a, b)) for b in range(16)] for a in range(16)], dtype=np.int64)
_hi, _lo = np.divmod(np.arange(N), 16)
# CHI[beta, alpha] = chi_beta(alpha)
CHI = 1 - 2 * (TRMUL[_hi[:, None], _hi[None, :]] ^ TRMUL[_lo[:, None], _lo[None, :]])
assert (CHI @ CHI == N * np.eye(N, dtype=np.int64)).all()


def walsh(f):
    return CHI @ np.asarray(f, dtype=np.int64)


def vec(x, y):
    return 16 * x + y


def line(m):
    """Nonzero points of the spread line with slope m (m = 16 means x = 0)."""
    if m == INF:
        return [vec(0, a) for a in range(1, 16)]
    return [vec(a, mul(m, a)) for a in range(1, 16)]


def annihilator_slope(m):
    """Slope of L_m^perp under the trace form: (x,y).(a,ma) = Tr(a(x+my)) = 0 for all a iff x = m y."""
    # L_m^perp = {(m y, y)}, which is L_{1/m} (and L_0^perp = L_inf, L_inf^perp = L_0).
    if m == 0:
        return INF
    if m == INF:
        return 0
    return inv(m)


def coset(a):
    """Index in Z/3 of the coset a C_5 of the order-5 subgroup C_5 = <t^3> of F^*."""
    return LOG[a] % 3


# ---------------------------------------------------------------- seeds
SEED_SLOPES = (0, 1, 2, 10)          # 2 = t, 10 = t^3 + t


def orbit_signing_to_z(slopes, negative):
    """The parametrization of the original notes: on each line (a, m a), order the three C_5-orbits
    by the smallest integer encoding the first coordinate, number the 12 orbits 0..11 in slope order,
    and put z = -1 on the orbits in `negative`."""
    z = np.zeros(N, dtype=np.int64)
    for k, m in enumerate(slopes):
        firsts = sorted(range(3), key=lambda j: min(a for a in range(1, 16) if coset(a) == j))
        for a in range(1, 16):
            idx = 3 * k + firsts.index(coset(a))
            z[line(m)[a - 1]] = -1 if idx in negative else 1
    return z


def is_seed(z):
    w = walsh(z)
    return set(np.unique(w).tolist()) <= {10, -6}


def seed_set(z):
    """B with H 1_B = 16 z + 96 delta_0."""
    w = walsh(z)
    assert set(np.unique(w).tolist()) == {10, -6}
    b = (w + 6) // 16
    assert ((w + 6) % 16 == 0).all() and b.sum() == 96
    hb = walsh(b)
    target = 16 * z.copy()
    target[0] += 96
    assert (hb == target).all()
    return b


BASE_SEED = orbit_signing_to_z(SEED_SLOPES, {0, 1, 3, 6, 11})


# ---------------------------------------------------------------- linear maps on V
def apply_matrix(Mx, v):
    """Mx = ((p, q), (r, s)) over F acts on column vectors: (a, b) -> (p a + q b, r a + s b)."""
    a, b = divmod(v, 16)
    (p, q), (r, s) = Mx
    return vec(mul(p, a) ^ mul(q, b), mul(r, a) ^ mul(s, b))


def transport(z, Mx):
    """z o M^{-1}: the seed moved by M on frequencies. Walsh values move to M^{-T} on vertices,
    so two-valuedness is preserved."""
    out = np.zeros(N, dtype=np.int64)
    for v in np.flatnonzero(z):
        out[apply_matrix(Mx, int(v))] = z[v]
    return out


def support_slopes(z):
    out = set()
    for m in list(range(16)) + [INF]:
        pts = line(m)
        if z[pts].any():
            assert (z[pts] != 0).all()
            out.add(m)
    return frozenset(out)


# ---------------------------------------------------------------- assembly
def assemble(seeds, translations=(0, 0, 0)):
    """Adjacency matrix (uint8, 768 x 768) for seeds z_0, z_1, z_2 on edges {0,1}, {1,2}, {2,0}.
    translations[i] = c replaces B_i by B_i + c (also allowed: H of the shifted set is 16 z chi_c + 96 delta)."""
    supports = [support_slopes(z) for z in seeds]
    assert all(len(s) == 4 for s in supports)
    assert not (supports[0] & supports[1] or supports[1] & supports[2] or supports[0] & supports[2])
    A = np.zeros((3 * N, 3 * N), dtype=np.uint8)
    xs = np.arange(N)
    xor = xs[:, None] ^ xs[None, :]
    for i in range(3):
        T = supports[i - 1] | supports[i]
        conn = np.zeros(N, dtype=np.uint8)
        for m in T:
            conn[line(annihilator_slope(m))] = 1
        assert conn.sum() == 120
        A[i * N:(i + 1) * N, i * N:(i + 1) * N] = conn[xor]
        b = seed_set(seeds[i])
        b = b[xs ^ translations[i]]
        j = (i + 1) % 3
        A[i * N:(i + 1) * N, j * N:(j + 1) * N] = b[xor]
        A[j * N:(j + 1) * N, i * N:(i + 1) * N] = b[xor].T
    return A


SHEARS = (0, 4, 13)


def shear_graph():
    """Three-orbit graph: the base seed sheared by (a, b) -> (a, b + h a), h = 0, 4, 13 (on frequencies)."""
    return assemble([transport(BASE_SEED, ((1, 0), (h, 1))) for h in SHEARS])


M_FREQ = ((1, 0), (5, 6))           # order 3 on frequencies


def cayley_seeds():
    z1 = transport(BASE_SEED, M_FREQ)
    z2 = transport(z1, M_FREQ)
    assert (transport(z2, M_FREQ) == BASE_SEED).all()
    return [BASE_SEED, z1, z2]


def cayley_graph():
    return assemble(cayley_seeds())


def vertex_map():
    """R = M^{-T} with respect to the trace form: Tr(<R v, M beta>) = Tr(<v, beta>)."""
    (p, q), (r, s) = M_FREQ
    det = mul(p, s) ^ mul(q, r)
    di = inv(det)
    # M^{-1} = det^{-1} ((s, q), (r, p)) in characteristic 2; transpose it.
    return ((mul(di, s), mul(di, r)), (mul(di, q), mul(di, p)))


def spectral_certificate(seeds):
    """Exact block-diagonalization: eigenvalues of A restricted to each character beta (3 x 3 block)."""
    supports = [support_slopes(z) for z in seeds]
    slope_of = {}
    for m in list(range(16)) + [INF]:
        for v in line(m):
            slope_of[v] = m
    eig = {}
    for beta in range(1, N):
        m = slope_of[beta]
        D = np.zeros((3, 3), dtype=np.int64)
        for i in range(3):
            D[i, i] = 8 if m in (supports[i - 1] | supports[i]) else -8
        for i in range(3):
            j = (i + 1) % 3
            D[i, j] = D[j, i] = 16 * seeds[i][beta]
        ev = np.round(np.linalg.eigvalsh(D.astype(float))).astype(int)
        assert np.allclose(np.linalg.eigvalsh(D.astype(float)), ev)
        for e in ev:
            eig[e] = eig.get(e, 0) + 1
    # beta = 0 block: diag 120, off-diagonal 96
    for e in (312, 24, 24):
        eig[e] = eig.get(e, 0) + 1
    return eig


if __name__ == '__main__':
    import sys
    sys.path.insert(0, __file__.rsplit('/', 1)[0])
    from verify_srg import check_srg
    print('base seed is two-valued:', is_seed(BASE_SEED), ' support slopes', sorted(support_slopes(BASE_SEED)))
    for name, seeds in (('shear', [transport(BASE_SEED, ((1, 0), (h, 1))) for h in SHEARS]),
                        ('cayley', cayley_seeds())):
        print(name, 'supports', [sorted(support_slopes(z)) for z in seeds])
        print(name, 'spectrum from 3x3 blocks', dict(sorted(spectral_certificate(seeds).items())))
        A = assemble(seeds)
        check_srg(A, (768, 312, 136, 120))
        check_srg(1 - A - np.eye(768, dtype=np.uint8), (768, 455, 262, 280))
        print(name, 'srg(768,312,136,120) and complement verified on all pairs')


# ---------------------------------------------------------------- the norm description
# N: F_16 -> F_4, N(a) = a^5 = a * a^4. F_4 = {0, 1, w, w^2} with w = t^2 + t (encoded 6), w^2 = 7.
OMEGA = 6
F4 = (0, 1, 6, 7)


def norm(a):
    return 0 if a == 0 else EXP[(5 * LOG[a]) % 15]


assert sorted({norm(a) for a in range(16)}) == list(F4) and mul(OMEGA, OMEGA) == 7 and mul(OMEGA, 7) == 1


def norm_seed(h, mu0, nu):
    """Seed on the slopes h + F_4: z(a, (h + mu) a) = eps_mu * (-1)^[N(a) = nu], eps_mu = -1 iff mu = mu0."""
    assert mu0 in F4 and nu in F4[1:]
    z = np.zeros(N, dtype=np.int64)
    for mu in F4:
        eps = -1 if mu == mu0 else 1
        for a in range(1, 16):
            z[vec(a, mul(h ^ mu, a))] = eps * (-1 if norm(a) == nu else 1)
    return z


COSET_REPS = (0, 2, 8, 10)         # 0, t, t^3, t^3+t represent the four additive cosets of F_4 = {0,1,6,7} in F_16


def parity_lemma():
    """For every w in F_16 and every 2-set {0, c} (c in F_4^*) scaled by any C in F_4^*, the number of m in F_4
    with N(w + m) in C*{0, c}... we check the statement used in the proof: for all w in F_16, n in F_4^*,
    #{m in F_4 : N(w + m) in {0, n}} is even; and for y != 0, N(x + m y) = N(y) N(x/y + m)."""
    for w in range(16):
        for n in F4[1:]:
            assert sum(norm(w ^ m) in (0, n) for m in F4) % 2 == 0
    for x in range(16):
        for y in range(1, 16):
            for m in F4:
                assert norm(x ^ mul(m, y)) == mul(norm(y), norm(mul(x, inv(y)) ^ m))
    return True


def gauss_period_lemma():
    """For u != 0: sum_{N(a) = nu} (-1)^Tr(a u) = -3 if N(u) nu = 1, else 1."""
    for nu in F4[1:]:
        for u in range(1, 16):
            s = sum((-1) ** tr(mul(a, u)) for a in range(1, 16) if norm(a) == nu)
            assert s == (-3 if mul(norm(u), nu) == 1 else 1), (nu, u, s)
    return True


def norm_threeorbit(cosets=(0, 2, 8), mu0s=(0, 0, 0), nus=(1, 1, 1), translation=0):
    return assemble([norm_seed(h, m, n) for h, m, n in zip(cosets, mu0s, nus)], translations=(translation, 0, 0))


def norm_cayley_seeds(c=2, mu0=0, nu=1):
    """M(a, b) = (a, c a + w b) has order 3 (c not in F_4) and moves the slope coset F_4 to c + F_4 and c w^2 + F_4."""
    M = ((1, 0), (c, OMEGA))
    z0 = norm_seed(0, mu0, nu)
    z1 = transport(z0, M)
    z2 = transport(z1, M)
    assert (transport(z2, M) == z0).all()
    return [z0, z1, z2], M
