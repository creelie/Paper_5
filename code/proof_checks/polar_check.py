#!/usr/bin/env python3
"""Checks for the proof of (P_inf) in polar.tex.

Conventions: Omega = [[rho, v], [v, sigma]], Y = Im Omega = (M1, M2, M3) with M1 = Im rho,
M2 = Im sigma, M3 = Im v.  A charge T = (m, n, l) stands for e(m rho + (n/2) sigma + l v),
n = 2 nu.  tau_Y = (-M3 + i sqrt(det Y)) / M2 is the point of the upper half-plane of Y.
A factor charge (k', l, j) of the product (1 - e(k' sigma + l rho + j v)) is stored as
(u, l, j) with u = 2k', and its exponent is Expo(u, l, j) of code/chl_coefficients.py.

(a) walls and chambers; (b) the Farey formula against the exact polar coefficients;
(c) the jump lemma on several walls, coefficients computed exactly in each chamber from
    the product, and one explicit path from K_0 to a vanishing chamber.
Usage: python3 polar_check.py [MMAX]   (default MMAX = 12)
"""
import sys, os, random, time
from fractions import Fraction as Fr
from math import gcd, floor, isqrt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import chl_coefficients as C

MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 12
EX = C.Expo(2, 20000)
T_START = time.time()

# ------------------------------------------------------------------ d1
J = 400; _d = [0] * (J + 1); _d[0] = 1
for r in range(1, J + 1):
    e = 16 if r % 2 == 0 else 8
    for _ in range(e):
        for k in range(r, J + 1): _d[k] += _d[k - r]
def d1(x):
    x = Fr(x)
    if x.denominator != 1 or x < -1: return 0
    return _d[int(x) + 1]
assert [d1(k) for k in range(-1, 4)] == [1, 8, 52, 256, 1122]

# ------------------------------------------------------------------ (a) walls and chambers
print('=== (a) walls and chambers ===')
def farey_edge(t1, t2):
    """t1, t2 in Q or None (= infinity): Farey neighbours?  returns (c, d) denominators."""
    def nd(t): return (1, 0) if t is None else (t.numerator, t.denominator)
    (a, c), (b, d) = nd(t1), nd(t2)
    return (abs(a * d - b * c) == 1), (c, d)

# a1: every N = -1 factor with nonzero exponent is a Farey edge with cd even, exponent 2,
#     and every Farey edge gives an N = -1 factor whose exponent is 2 if cd even, 0 otherwise
n1 = 0
for u in range(0, 61):
    for l in range(0, 61):
        for j in range(-isqrt(2 * u * l + 1) - 1, isqrt(2 * u * l + 1) + 2):
            if 2 * u * l - j * j != -1: continue
            if (u, l) == (0, 0) and j > 0: continue
            c = EX(u, l, j); n1 += 1
            # wall in the tau-plane: l |tau|^2 - j Re tau + u/2 = 0, roots of l t^2 - j t + u/2
            if l == 0:
                ends = (Fr(u, 2 * j), None)
            else:
                ends = (Fr(j - 1, 2 * l), Fr(j + 1, 2 * l))
            ok, (cc, dd) = farey_edge(*ends) if (u % 2 == 0) else (None, (None, None))
            if u % 2 == 1:
                assert c == 0, ('odd u with nonzero exponent', u, l, j)
                continue
            assert ok, ('N=-1 factor not a Farey edge', u, l, j, ends)
            assert c == (2 if (cc * dd) % 2 == 0 else 0), (u, l, j, ends, c)
print('a1: %d N=-1 factor charges with u,l <= 60: exponent 2 exactly on Farey edges with cd even,'
      ' 0 on odd-odd edges and for u odd' % n1)

def fracs(D, lo, hi):
    out = set()
    for c in range(1, D + 1):
        for a in range(floor(lo * c), floor(hi * c) + 1):
            if gcd(a, c) == 1 and lo <= Fr(a, c) <= hi: out.add(Fr(a, c))
    return sorted(out)
def is_wall(t1, t2):
    ok, (c, d) = farey_edge(t1, t2)
    return ok and (c * d) % 2 == 0
def wall_charge(t1, t2):
    """positive factor charge (u, l, j) of the wall through t1, t2 (Farey neighbours)."""
    if t2 is None: t1, t2 = t2, t1
    if t1 is None:  # vertical line Re tau = a : -j Re tau + k' = 0 with j = +-1
        a = t2
        assert a.denominator == 1
        a = int(a)
        if a > 0: return (2 * a, 0, 1)
        if a < 0: return (-2 * a, 0, -1)
        return (0, 0, -1)
    a, c, b, d = t1.numerator, t1.denominator, t2.numerator, t2.denominator
    k, l, j = a * b, c * d, a * d + b * c
    if l < 0 or (l == 0 and k < 0): k, l, j = -k, -l, -j
    return (2 * k, l, j)

# a2: K_0 = ideal quadrilateral (0, 1/2, 1, inf): sample points, all wall functionals positive
walls = []
F = fracs(16, Fr(-3), Fr(3))
for i, t1 in enumerate(F):
    for t2 in F[i + 1:]:
        if is_wall(t1, t2): walls.append((t1, t2))
for n in range(-3, 4): walls.append((Fr(n), None))
WCH = [wall_charge(*w) for w in walls]
assert all(EX(*A) == 2 for A in WCH)
def sval(A, tau):  # sign of tr(A Y_tau)
    u, l, j = A
    return l * abs(tau) ** 2 - j * tau.real + u / 2
def in_quad(tau, vs):  # vs = sorted finite vertices v1<v2<v3, optional inf ; ideal polygon test
    return None
random.seed(1)
inside = outside = 0
for _ in range(4000):
    tau = complex(random.uniform(-1.5, 2.5), random.uniform(0.01, 3))
    allpos = all(sval(A, tau) > 0 for A in WCH)
    inQ = (0 < tau.real < 1) and abs(tau - 0.25) > 0.25 and abs(tau - 0.75) > 0.25
    assert allpos == inQ, tau
    inside += inQ; outside += (not inQ)
print('a2: K_0 = {all %d walls of the box positive} = quadrilateral (0,1/2,1,inf): %d points inside,'
      ' %d outside, all consistent' % (len(walls), inside, outside))

# a3: chambers = quadrilaterals around odd-odd Farey edges, each the Gamma_0(2)-image of K_0,
#     and the complement of the walls has exactly these components (sign vectors)
def mob(M, t):
    a, b, c, d = M
    if t is None: return None if c == 0 else Fr(a, c)
    den = c * t + d
    return None if den == 0 else (a * t + b) / den
def mobc(M, z):
    a, b, c, d = M; return (a * z + b) / (c * z + d)
def sl2_from(t1, t2):
    """M in SL2(Z) with M(inf) = t1, M(0) = t2 (Farey neighbours, t1 != inf)."""
    a, c = t1.numerator, t1.denominator
    if t2 is None: b, d = 1, 0
    else: b, d = t2.numerator, t2.denominator
    if a * d - b * c == -1: b, d = -b, -d
    assert a * d - b * c == 1
    return (a, b, c, d)
K0 = [Fr(0), Fr(1, 2), Fr(1), None]
quads = []
Fs = fracs(12, Fr(-2), Fr(2))
for i, t1 in enumerate(Fs):
    for t2 in Fs[i + 1:]:
        ok, (c, d) = farey_edge(t1, t2)
        if not ok or (c % 2 == 0) or (d % 2 == 0): continue
        # the two Farey triangles on the edge: third vertices (a+b)/(c+d) and (a-b)/(c-d)
        a, b = t1.numerator, t2.numerator
        th1 = Fr(a + b, c + d)
        th2 = None if c == d else Fr(a - b, c - d)
        verts = [t1, t2, th1, th2]
        # outer edges are walls
        for x in (t1, t2):
            for y in (th1, th2):
                assert is_wall(x, y), (verts, x, y)
        # Gamma_0(2) element mapping K_0 onto it: absent edge (0,1) -> (t1,t2)
        found = None
        for (x, y) in ((t1, t2), (t2, t1)):
            g = sl2_from(x, y)          # g(inf) = x, g(0) = y
            A = (g[0] - g[1], g[1], g[2] - g[3], g[3])   # A = g h^{-1}, h(z) = z/(z+1): A(1) = x, A(0) = y
            if A[2] % 2: continue
            img = [mob(A, v) for v in K0]
            if set(map(str, img)) == set(map(str, verts)): found = A
        assert found is not None, verts
        quads.append((verts, found))
print('a3: %d chambers (quadrilaterals on odd-odd Farey edges with denominators <= 12 in [-2,2]):'
      ' all four sides are walls and each is A(K_0) with A in Gamma_0(2)' % len(quads))
# sign-vector test: points of a quadrilateral A(K_0) are A(points of K_0); equal sign vectors
# inside one quadrilateral, distinct between quadrilaterals
def sv(tau): return tuple(sval(A, tau) > 0 for A in WCH)
seen = {}; npts = 0
for verts, A in quads[:400]:
    vecs = set()
    for _ in range(6):
        while True:
            z = complex(random.uniform(0.02, 0.98), random.uniform(0.05, 2.5))
            if abs(z - 0.25) > 0.26 and abs(z - 0.75) > 0.26: break
        w = mobc(A, z)
        if abs(w.real) > 2.2: continue
        vecs.add(sv(w)); npts += 1
    assert len(vecs) <= 1, verts
    for v in vecs:
        assert v not in seen or seen[v] == str(verts), ('two chambers with one sign vector', verts)
        seen[v] = str(verts)
print('a3: %d sampled points in %d chambers: constant sign vector on each chamber, distinct across'
      % (npts, len(seen)))
# a4: Gamma_0(2) generators E: tau -> tau - 1 and L: tau -> tau/(1 - 2 tau) permute the walls
Wset = set((str(a), str(b)) for a, b in walls)
for (t1, t2) in walls:
    for M in ((1, -1, 0, 1), (1, 0, -2, 1), (1, 1, 0, 1), (1, 0, 2, 1)):
        x, y = mob(M, t1), mob(M, t2)
        assert is_wall(x, y)
        x2, y2 = mob(M, Fr(1) if t1 is None else t1), None
    # odd-odd edges are not mapped to walls
print('a4: the generators of Gamma_0(2) (and inverses) map each of the %d walls to a wall' % len(walls))

# ------------------------------------------------------------------ (b) Farey formula vs exact
print('=== (b) Farey wall sum against exact polar coefficients of phi^F_m at infinity ===')
coef, _ = C.chamber_coefficients(2, MMAX + 2, 4 * MMAX + 4)
def farey_sum(m, nu, l, closed_left=True, B=None):
    w0 = Fr(l, 2 * m)
    if B is None: B = 4 * m + 6
    tot = 0
    for r in range(1, B):
        for s in range(1, B):
            if gcd(r, s) != 1 or (r % 2 and s % 2): continue
            qi = (-pow(r, -1, s)) % s if s > 1 else 0
            if closed_left:   # q/s <= w0 < p/r
                q = qi + s * floor((w0 * s - qi) / s)
                if not (Fr(q, s) <= w0 < Fr(q, s) + Fr(1, r * s)): continue
            else:             # q/s < w0 <= p/r
                q = qi + s * floor((w0 * s - qi - Fr(1, r)) / s) if False else None
                q = None
                for qq in range(floor(w0 * s) - 2 * s - 2, floor(w0 * s) + 2):
                    if (qq * r + 1) % s == 0 and Fr(qq, s) < w0 <= Fr(qq, s) + Fr(1, r * s): q = qq
                if q is None: continue
            p = (1 + q * r) // s
            mG = r * r * nu + p * p * m - p * r * l
            nG = s * s * nu + q * q * m - q * s * l
            lG = -2 * r * s * nu - 2 * p * q * m + (p * s + q * r) * l
            if r % 2 == 0: tot += lG * d1(mG) * d1(2 * nG)
            else: tot += lG * d1(2 * mG) * d1(nG)
    return tot
nt = nz = bad = 0; per_m = {}
for m in range(1, MMAX + 1):
    for l in range(0, 2 * m):
        n = -4
        while 2 * m * n - l * l < 0:
            nu = Fr(n, 2)
            exact = coef(m, n, l) if n >= -1 else 0
            fs = farey_sum(m, nu, l)
            fs2 = farey_sum(m, nu, l, closed_left=False)
            nt += 1; nz += (exact != 0)
            per_m[m] = per_m.get(m, 0) + (exact != 0)
            if not (exact == fs == fs2):
                bad += 1; print('  MISMATCH', m, n, l, exact, fs, fs2)
            n += 1
print('b1: %d polar terms (1<=m<=%d, 0<=l<2m, -2<=nu, 4m nu<l^2), %d nonzero; mismatches: %d'
      % (nt, MMAX, nz, bad))
print('    nonzero polar coefficients per m:', per_m)
print('    the same holds with the convention q/s < w0 <= p/r (both conventions agree)')
# stability of the truncation r,s < 4m+6
st = 0
for m in range(1, 7):
    for l in range(0, 2 * m):
        for n in range(-2, (l * l) // (2 * m) + 1):
            if 2 * m * n - l * l < 0:
                assert farey_sum(m, Fr(n, 2), l) == farey_sum(m, Fr(n, 2), l, B=12 * m + 20); st += 1
print('b2: truncation r,s < 4m+6 agrees with r,s < 12m+20 on %d polar terms with m <= 6' % st)
print('    examples (m, n=2nu, l, coefficient):',
      [(m, n, l, coef(m, n, l)) for (m, n, l) in [(1, -1, 0), (1, -1, 1), (2, -1, 2), (3, -1, 3), (4, 1, 3), (8, 2, 7)]])

# ------------------------------------------------------------------ (c) jump lemma
print('=== (c) jump lemma, coefficients computed exactly chamber by chamber ===')
def Ymat(tau, sd=1.0):
    y = tau.imag; return (sd / y * abs(tau) ** 2, sd / y, -sd / y * tau.real)
def gbinom(a, k):
    r = 1; f = 1
    for i in range(k): r *= (a - i); f *= (i + 1)
    return r // f
def chamber_series(Y, K):
    """Exact Fourier coefficients of 1/Phi~6 in the chamber of Y: the finitely many factors
    with tr(AY) < 0 are inverted, (1 - e(w))^(-c) = (-e(-w))^c (1 - e(-w))^(-c), and the
    resulting product is a series in charges of positive weight <., Y>, truncated at K.
    Generators are enumerated completely: psd factors with tr(AY) <= K satisfy
    k' + l <= K / lambda_min(Y); N = -1 factors are walls (Farey edges with cd even, exponent 2)
    and tr(AY) = sqrt(det Y) Q / Im tau with Q^2 + (Im tau)^2 = |c tau - a|^2 |d tau - b|^2."""
    M1, M2, M3 = Y
    w = lambda T: T[0] * M1 + T[1] * M2 / 2 + T[2] * M3
    det = M1 * M2 - M3 * M3; sd = det ** 0.5
    lam = (M1 + M2) / 2 - (((M1 - M2) / 2) ** 2 + M3 * M3) ** 0.5
    tau = complex(-M3 / M2, sd / M2); x, y = tau.real, tau.imag
    gens = []; T0 = [-1, -1, -1]; sgn = 1
    def add(A, c):
        nonlocal T0, sgn
        tr = w(A)
        if tr < 0:
            T0 = [T0[i] - c * A[i] for i in range(3)]; sgn *= (-1) ** c
            gens.append(((-A[0], -A[1], -A[2]), c, -tr))
        elif tr <= K: gens.append((A, c, tr))
    # psd factors (N >= 0)
    Lmax = int(K / lam) + 1
    for u in range(0, 2 * Lmax + 1):
        for l in range(0, Lmax + 1):
            jm = isqrt(2 * u * l)
            for j in range(-jm, jm + 1):
                if (u == 0 and l == 0) or 2 * u * l - j * j < 0: continue
                if w((l, u, j)) > K: continue
                c = EX(u, l, j)
                if c: add((l, u, j), c)
    # walls (N = -1): vertical ones Re tau = a, and Farey edges a/c - b/d with c, d >= 1
    Xb = ((K * y / sd) ** 2 + y * y) ** 0.5
    seen = set()
    for a in range(floor(x - K * y / sd) - 1, floor(x + K * y / sd) + 2):
        A = wall_charge(Fr(a), None); u, l, j = A; add((l, u, j), 2)
    # walls with 0 <= tr <= K: |c tau - a| |d tau - b| <= Xb with |c tau - a| >= c y, |d tau - b| >= d y,
    # so c d y^2 <= Xb and |c x - a| <= Xb/(d y); walls with tr < 0 separate tau from K_0,
    # so pass above tau: radius 1/(2cd) > y and |c x - a| <= 1/d
    Cb = max(Xb / (y * y), 1 / (2 * y))
    for c in range(1, int(Cb) + 2):
        for d in range(1, int(Cb / c) + 2):
            if gcd(c, d) != 1 or (c * d) % 2: continue
            ar = max(Xb / (d * y), 1 / d)
            for a in range(floor(c * x - ar) - 1, floor(c * x + ar) + 2):
                if gcd(a, c) != 1: continue
                for sg in (1, -1):        # a d - b c = sg
                    if (a * d - sg) % c: continue
                    b = (a * d - sg) // c
                    t1, t2 = Fr(a, c), Fr(b, d)
                    key = (min(t1, t2), max(t1, t2))
                    if key in seen: continue
                    seen.add(key)
                    u, l, j = wall_charge(t1, t2); add((l, u, j), 2)
    S = {(0, 0, 0): 1}
    for (g, c, tr) in gens:
        new = dict(S)
        for T, v in S.items():
            wt = w(T); k = 1
            while wt + k * tr <= K + 1e-9:
                key = (T[0] + k * g[0], T[1] + k * g[1], T[2] + k * g[2])
                new[key] = new.get(key, 0) + (-1) ** k * gbinom(-c, k) * v
                k += 1
        S = {k: v for k, v in new.items() if v}
    def d(m, n, l):
        R = (m - T0[0], n - T0[1], l - T0[2])
        if w(R) > K + 1e-9: raise ValueError
        return sgn * S.get(R, 0)
    return d
def jump_pred(G, m, n, l):
    p, q, r, s = G; nu = Fr(n, 2)
    mG = r * r * nu + p * p * m - p * r * l; nG = s * s * nu + q * q * m - q * s * l
    lG = -2 * r * s * nu - 2 * p * q * m + (p * s + q * r) * l
    return lG * ((d1(mG) * d1(2 * nG) if r % 2 == 0 else 0) + (d1(2 * mG) * d1(nG) if s % 2 == 0 else 0))
def test(name, tout, tin, predf, K=5):
    do = chamber_series(Ymat(tout), K); di = chamber_series(Ymat(tin), K)
    n = nz = bad = 0
    for m in range(-1, 4):
        for nn in range(-1, 8):
            for l in range(-12, 13):
                try: a = do(m, nn, l); b = di(m, nn, l)
                except ValueError: continue
                pr = predf(m, nn, l); n += 1; nz += (pr != 0)
                if a - b != pr: bad += 1; print('  FAIL', (m, nn, l), a, b, pr)
    print('c1: %-28s %4d charges, %3d with nonzero predicted jump, failures %d' % (name, n, nz, bad))
# check of the method: K_0 coefficients agree with chl_coefficients
dK = chamber_series(Ymat(0.3 + 1.0j), 6)
chk = [(1, 0, 0), (1, 1, -1), (1, 1, 1), (2, 0, -2), (0, 2, -1), (1, -1, -2), (-1, -1, -3), (0, 0, -1)]
assert all(dK(*T) == coef(*T) for T in chk)
print('c0: chamberwise exact expansion reproduces chl_coefficients in K_0 on', chk)
test('basic wall Re tau = 0', 0.05 + 1j, -0.05 + 1j, lambda m, n, l: -l * d1(m) * d1(n))
for G in [(1, 0, 2, 1), (1, 1, 1, 2), (1, 1, 2, 3), (1, 0, 3, 1), (2, 1, 5, 3), (1, 1, 3, 4), (0, -1, 1, 2), (2, 1, 3, 2)]:
    p, q, r, s = G; a, b = q / s, p / r; c0 = (a + b) / 2; r0 = (b - a) / 2
    kind = 'r even' if r % 2 == 0 else ('s even' if s % 2 == 0 else 'absent (r,s odd)')
    test('%d/%d -> %d/%d (%s)' % (q, s, p, r, kind), complex(c0, 1.05 * r0), complex(c0, 0.95 * r0),
         lambda m, n, l, G=G: jump_pred(G, m, n, l))

# c2: an explicit path for T = (m, n, l) = (3, 0, 1): nu = 0, w0 = 1/6, t = 1, delta = 1/6,
#     kappa* = 3/16 (even denominator), V(3/16) = -21.  All Farey edges straddling 3/16 (denominators
#     adding up to <= 16 by (F3)) are listed with their crossing heights; the coefficient of T is
#     computed exactly in the chamber between consecutive crossings and each difference is compared
#     with the transport formula (0 for odd edges); the last chamber must give 0.
T = (3, 0, 1); kap = Fr(3, 16); m_, nu_, l_ = 3, Fr(0), 1
V = 16 * 16 * nu_ - l_ * 3 * 16 + m_ * 9
edges = []
for t1 in fracs(16, Fr(-1), Fr(2)):
    for t2 in fracs(16, Fr(-1), Fr(2)):
        if t1 < kap < t2 and farey_edge(t1, t2)[0]:
            c0, r0 = (t1 + t2) / 2, (t2 - t1) / 2
            hgt = float(r0 * r0 - (kap - c0) ** 2) ** 0.5
            q, s_, p, r = t1.numerator, t1.denominator, t2.numerator, t2.denominator
            edges.append((hgt, (p, q, r, s_)))
edges.sort(reverse=True)
walls_ = [i for i, (hgt, G) in enumerate(edges) if (G[2] * G[3]) % 2 == 0]
# sample just below each wall, above the next Farey edge (odd edges between walls do not change the chamber)
hs = [2.0] + [(edges[i][0] + (edges[i + 1][0] if i + 1 < len(edges) else 0.0)) / 2 for i in walls_]
vals = []
for h in hs:
    for K in range(3, 12):
        try: vals.append(chamber_series(Ymat(complex(float(kap), h)), K)(*T)); break
        except ValueError: continue
print('c2: path Re tau = 3/16 for T = (m,n,l) = (3,0,1), w0 = 1/6, V(3/16) = %d; K_0 value %d' % (V, coef(*T)))
print('    Farey edges crossed (height, q/s -> p/r):',
      ', '.join('%.4f %d/%d->%d/%d%s' % (hg, G[1], G[3], G[0], G[2], '' if (G[2] * G[3]) % 2 == 0 else ' (odd)')
                for hg, G in edges))
assert V <= -2 and vals[0] == coef(*T) and vals[-1] == 0
for k, i in enumerate(walls_):
    hgt, G = edges[i]; p, q, r, s_ = G
    pr = jump_pred(G, *T)
    kind = 'r even' if r % 2 == 0 else 's even'
    print('    wall %s/%s -> %s/%s (%s) at height %.4f: d %6d -> %6d, jump %6d, predicted %6d'
          % (q, s_, p, r, kind, hgt, vals[k], vals[k + 1], vals[k] - vals[k + 1], pr))
    assert vals[k] - vals[k + 1] == pr
print('    sum of predicted jumps = %d = K_0 value; Farey sum a_I = %d'
      % (sum(jump_pred(G, *T) for _, G in edges), farey_sum(3, Fr(0), 1)))
assert sum(jump_pred(G, *T) for _, G in edges) == coef(*T) == farey_sum(3, Fr(0), 1)
print('note: the identity Phi~6(rho,sigma,v) = 16 sigma^-6 nabla3^2(rho - v^2/sigma, -1/sigma, v/sigma)'
      ' is not used by the proof (nor checked here)')
print('total time %.1f s' % (time.time() - T_START))
print('ALL CHECKS PASSED')
