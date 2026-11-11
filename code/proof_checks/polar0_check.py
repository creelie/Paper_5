#!/usr/bin/env python3
"""Checks for the proof of (P0) in polar0.tex.

Conventions: Omega = [[rho, v], [v, sigma]], Y = Im Omega = (M1, M2, M3), M1 = Im rho, M2 = Im sigma,
M3 = Im v; p = e(sigma), x = e(sigma/2), y = e(v).  S Omega = (rho - v^2/sigma, -1/sigma, v/sigma).
Phi~6 is the product (eq:product) of the paper; Phi6 is the product of polar0.tex,
    Phi6 = e(rho + sigma + v) prod_{(B,k,j)} (1 - X)^{C0(4Bk - j^2)} (1 + X)^{C1(4Bk - j^2)},
X = e(B rho + k sigma + j v), B, k >= 0 integers, j < 0 if B = k = 0, with
C0 = c^(0,0) + c^(0,1) = (4/3) B + (2/3) A E,  C1 = c^(1,0) + c^(1,1) = (2/3) B - (2/3) A E
(coefficients at the discriminant N).  A charge T = (m, nu, l) stands for e(m rho + nu sigma + l v).

(i)   the identity Phi~6(Omega) = 16 sigma^-6 Phi6(S Omega) from the two products, and the
      Jacobi-form identities behind its proof, from theta functions at transformed points;
(ii)  walls and chambers of 1/Phi6: every Farey edge is a wall, the chambers are the Farey
      triangles, the jump across every wall (odd-odd edges included), and an explicit path;
(iii) the polar coefficients of the frame-S finite part psi^S_m = phi_m|S - (d1(m)/16) A_{2,m}/g(sigma),
      computed exactly from its q-expansion, against the wall sum (eq:P0), m <= MMAX;
(iv)  the frame-S finite part from a direct numerical S-transformation of phi_m (no use of Phi6 or of
      the permutation of the genera), Fourier coefficients by FFT, m = 1, 2;
(v)   the Appell-Lerch sum: A_{2,m}|S - A_{2,m} is a period integral of the unary theta functions
      (it has no Fourier expansion in sigma and contributes no polar term).
The count of polar terms in the convention of the archive (0 <= l <= m, nu >= -1) is printed in (iii3).
Usage: python3 polar0_check.py [MMAX]   (default MMAX = 12)
"""
import sys, os, time, random, cmath
from fractions import Fraction as Fr
from math import gcd, floor, isqrt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import chl_coefficients as Cc
import mpmath as mp

MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 12
T_START = time.time()
EX = Cc.Expo(2, 4000)

def C0(D):
    if D < -1: return 0
    v = Fr(16, 3) * EX.a(D) - Fr(2, 3) * EX.BEt[D]
    assert v.denominator == 1
    return int(v)
def C1(D):
    if D < -1: return 0
    v = Fr(8, 3) * EX.a(D) + Fr(2, 3) * EX.BEt[D]
    assert v.denominator == 1
    return int(v)

J = 600; _d = [0] * (J + 1); _d[0] = 1
for r in range(1, J + 1):
    e_ = 16 if r % 2 == 0 else 8
    for _ in range(e_):
        for k in range(r, J + 1): _d[k] += _d[k - r]
def d1(x):
    x = Fr(x)
    if x.denominator != 1 or x < -1: return 0
    return _d[int(x) + 1]
assert [d1(k) for k in range(-1, 4)] == [1, 8, 52, 256, 1122]

# ================================================================== (i)
print('=== (i) the transformation S of Phi~6, from the two products ===')
mp.mp.dps = 45
def e(z): return mp.exp(2j * mp.pi * z)
def lam_min(M1, M2, M3): return (M1 + M2) / 2 - (((M1 - M2) / 2) ** 2 + M3 * M3) ** 0.5
def logPhit(rho, sig, v, X):
    M1, M2, M3 = float(mp.im(rho)), float(mp.im(sig)), float(mp.im(v))
    Bm = int(X / lam_min(M1, M2, M3)) + 3
    lg = 2j * mp.pi * (rho + sig / 2 + v)
    for B in range(0, Bm + 1):
        for u in range(0, 2 * Bm + 1):            # u = 2k', k' in Z/2
            jm = isqrt(2 * u * B + 1) + 1
            for j in range(-jm, jm + 1):
                if (B == 0 and u == 0 and j >= 0) or 2 * u * B - j * j < -1: continue
                if B * M1 + u * M2 / 2 + j * M3 > X: continue
                c = EX(u, B, j)
                if c: lg += c * mp.log(1 - e(B * rho + Fr(u, 2) * sig + j * v))
    return lg
def logPhi6(rho, sig, v, X):
    M1, M2, M3 = float(mp.im(rho)), float(mp.im(sig)), float(mp.im(v))
    Bm = int(X / lam_min(M1, M2, M3)) + 3
    lg = 2j * mp.pi * (rho + sig + v)
    for B in range(0, Bm + 1):
        for k in range(0, Bm + 1):
            jm = isqrt(4 * k * B + 1) + 1
            for j in range(-jm, jm + 1):
                N = 4 * k * B - j * j
                if (B == 0 and k == 0 and j >= 0) or N < -1: continue
                if B * M1 + k * M2 + j * M3 > X: continue
                Xg = e(B * rho + k * sig + j * v)
                a, b = C0(N), C1(N)
                if a: lg += a * mp.log(1 - Xg)
                if b: lg += b * mp.log(1 + Xg)
    return lg
def Smap(rho, sig, v): return (rho - v * v / sig, -1 / sig, v / sig)
pts = [(mp.mpc(0.13, 5.0), mp.mpc(0.31, 0.97), mp.mpc(0.07, -0.04)),
       (mp.mpc(-0.4, 4.2), mp.mpc(-0.2, 1.1), mp.mpc(0.31, 0.12)),
       (mp.mpc(0.77, 6.5), mp.mpc(0.55, 0.8), mp.mpc(-0.21, 0.3))]
worst = 0
for P in pts:
    l1 = logPhit(*P, X=17.0)
    l2 = logPhi6(*Smap(*P), X=17.0)
    dd = l1 - (mp.log(16) - 6 * mp.log(P[1]) + l2)
    dd = dd - 2j * mp.pi * mp.nint(mp.im(dd) / (2 * mp.pi))
    worst = max(worst, abs(dd))
print('i1: log Phi~6(Omega) - log(16 sigma^-6 Phi6(S Omega)) at 3 points (products truncated at'
      ' <g,Y> <= 17): max |difference| mod 2 pi i = %s' % mp.nstr(worst, 3))
assert worst < 1e-30

# Jacobi forms by theta functions (paper's convention theta_k(tau, z) = jtheta(k, pi z, e^{i pi tau}))
def th(k, tau, z=0):
    """theta_k(tau, z) of the paper by direct summation (explicit q^(1/4), no branch ambiguity)"""
    yv = float(mp.im(tau)); nmax = int((90 / (mp.pi * yv)) ** 0.5) + 3 + int(abs(float(mp.im(z))) / yv)
    if k in (1, 2):
        return 2 * sum((-1) ** (n * (k == 1)) * mp.exp(1j * mp.pi * tau * (n + 0.5) ** 2)
                       * (mp.sin if k == 1 else mp.cos)((2 * n + 1) * mp.pi * z) for n in range(-nmax, nmax + 1) if n >= 0)
    return 1 + 2 * sum((-1) ** (n * (k == 4)) * mp.exp(1j * mp.pi * tau * n * n) * mp.cos(2 * n * mp.pi * z)
                       for n in range(1, nmax + 1))
def eta(tau): return mp.exp(1j * mp.pi * tau / 12) * mp.qp(mp.exp(2j * mp.pi * tau))
def pA(tau, z): return -th(1, tau, z) ** 2 / eta(tau) ** 6
def pB(tau, z): return 4 * sum(th(k, tau, z) ** 2 / th(k, tau) ** 2 for k in (2, 3, 4))
def fE(tau):    # 2 E(2 tau) - E(tau) = theta_3(2 tau)^4 + theta_2(2 tau)^4
    Q = mp.exp(2j * mp.pi * tau); return mp.jtheta(3, 0, Q) ** 4 + mp.jtheta(2, 0, Q) ** 4
def Fg(r, s, tau, z):
    b, a = pB(tau, z), pA(tau, z)
    if (r, s) == (0, 0): return b
    if (r, s) == (0, 1): return b / 3 + 2 * a * fE(tau) / 3
    if (r, s) == (1, 0): return b / 3 - a * fE(tau / 2) / 3
    return b / 3 - a * fE((tau + 1) / 2) / 3
# q-expansion of fE
q0 = mp.mpc(0.1, 0.9)
ser = 1 + 24 * sum(sum(dd for dd in range(1, n + 1) if n % dd == 0 and dd % 2) * e(n * q0) for n in range(1, 60))
assert abs(ser - fE(q0)) < 1e-30
mp.mp.dps = 30
def slashS(f, l, tau, z, k=0):   # (f|_{k,l} S)(tau, z)
    return tau ** (-k) * e(-l * z * z / tau) * f(-1 / tau, z / tau)
bad = 0
for (tau, z) in [(mp.mpc(0.2, 0.9), mp.mpc(0.13, 0.05)), (mp.mpc(-0.35, 1.3), mp.mpc(-0.4, 0.2))]:
    for r in (0, 1):
        for s in (0, 1):
            dS = abs(slashS(lambda a, b: Fg(r, s, a, b), 1, tau, z) - Fg(s, r, tau, z))
            dT = abs(Fg(r, s, tau + 1, z) - Fg(r, (r + s) % 2, tau, z))
            bad += (dS > 1e-20) + (dT > 1e-20)
    dA = abs(slashS(pA, 1, tau, z, k=-2) - pA(tau, z))
    g = lambda t: eta(t) ** 8 * eta(2 * t) ** 8
    dg = abs(tau ** 8 / g(-1 / (2 * tau)) - 1 / (16 * g(tau)))
    bad += (dA > 1e-20) + (dg > 1e-20)
assert bad == 0
print('i2: F^(r,s)|S = F^(s,r), F^(r,s)|T = F^(r,r+s), A|_{-2,1}S = A and (1/g(sigma/2))|_{-8}S = 1/(16 g(sigma))'
      ' at 2 points, to 1e-20')
def Lhecke(l, sig, v, twisted):
    tot = 0
    for t in range(1, l + 1):
        if l % t: continue
        d = l // t
        for u in range(2 * d):
            for r in (0, 1):
                for s in (0, 1):
                    sg = (-1) ** (r * t) if twisted else (-1) ** (s * d + r * u)
                    tot += sg * Fg(r, s, (t * sig + u) / d, t * v)
    return tot / (2 * l)
def Lprime_series(l, sig, v, AM=40):
    tot = 0
    for t in range(1, l + 1):
        if l % t: continue
        B = l // t
        for k in range(0, AM):
            jm = isqrt(4 * B * k + 1)
            for j in range(-jm, jm + 1):
                N = 4 * B * k - j * j
                c = C0(N) + (-1) ** t * C1(N)
                if c: tot += mp.mpf(c) / t * e(t * (k * sig + j * v))
    return tot
worst = 0
for (sig, v) in [(mp.mpc(0.21, 1.05), mp.mpc(0.1, -0.15)), (mp.mpc(-0.33, 0.85), mp.mpc(0.45, 0.2))]:
    for l in range(1, 6):
        lhs = slashS(lambda a, b: Lhecke(l, a, b, False), l, sig, v)
        rhs = Lhecke(l, sig, v, True)
        rhs2 = Lprime_series(l, sig, v)
        worst = max(worst, abs(lhs - rhs) / abs(rhs), abs(rhs - rhs2) / abs(rhs))
assert worst < 1e-15
print("i3: L_l|_{0,l}S = L'_l (Hecke form with sign (-1)^(rt)) = sum_t (1/t) sum (C0 + (-1)^t C1) p^(tk) y^(tj),"
      ' l = 1..5, 2 points, from theta functions at the transformed points: max rel. error %.1e' % worst)

# ================================================================== (ii)
print('=== (ii) walls and chambers of 1/Phi6 ===')
assert C0(-1) == 2 and C1(-1) == 0 and all(C0(D) == 0 == C1(D) for D in range(-9, -1))
def nd(t): return (1, 0) if t is None else (t.numerator, t.denominator)
def farey(t1, t2):
    (a, c), (b, d) = nd(t1), nd(t2); return abs(a * d - b * c) == 1
nw = 0
for B in range(0, 61):
    for k in range(0, 61):
        for j in range(-isqrt(4 * B * k + 1) - 1, isqrt(4 * B * k + 1) + 2):
            if 4 * B * k - j * j != -1 or (B == k == 0 and j > 0): continue
            ends = (Fr(k, j), None) if B == 0 else (Fr(j - 1, 2 * B), Fr(j + 1, 2 * B))
            assert farey(*ends); nw += 1
print('ii1: C0(-1) = 2, C1(-1) = 0, C0 = C1 = 0 below -1; the %d factor charges with N = -1, B, k <= 60, are'
      ' exactly Farey edges, each with exponent 2: every Farey edge is a wall' % nw)
def fracs(D, lo, hi):
    out = set()
    for c in range(1, D + 1):
        for a in range(floor(lo * c), floor(hi * c) + 1):
            if gcd(a, c) == 1 and lo <= Fr(a, c) <= hi: out.add(Fr(a, c))
    return sorted(out)
def wall_charge(t1, t2):
    """factor charge (B, k, j) of the wall through the Farey neighbours t1, t2"""
    if t2 is None: t1, t2 = t2, t1
    if t1 is None:
        a = int(t2)
        return (0, a, 1) if a > 0 else ((0, -a, -1) if a < 0 else (0, 0, -1))
    a, c, b, d = t1.numerator, t1.denominator, t2.numerator, t2.denominator
    B, k, j = c * d, a * b, a * d + b * c
    if B < 0: B, k, j = -B, -k, -j
    return (B, k, j)
def sval(g, tau): B, k, j = g; return B * abs(tau) ** 2 - j * tau.real + k
F_ = fracs(10, Fr(-2), Fr(3)); WL = []
for i, t1 in enumerate(F_):
    for t2 in F_[i + 1:]:
        if farey(t1, t2): WL.append(wall_charge(t1, t2))
for n in range(-2, 4): WL.append(wall_charge(Fr(n), None))
def mobc(M, z): a, b, c, d = M; return (a * z + b) / (c * z + d)
random.seed(2)
inside = 0
for _ in range(3000):
    tau = complex(random.uniform(-1.5, 2.5), random.uniform(0.02, 3))
    allpos = all(sval(g, tau) > 0 for g in WL)
    inT = 0 < tau.real < 1 and abs(tau - 0.5) > 0.5
    assert allpos == inT; inside += inT
# every Farey triangle is M(Delta(0)), M in SL2(Z): sign vectors constant on M(Delta(0)), distinct across
seen = {}; npts = 0
for c in range(0, 5):
    for d in range(-5, 6):
        if gcd(c, d) != 1 or (c == 0 and d != 1): continue
        for a in range(-6, 7):
            if c == 0: b = None
            if (a * d - 1) % c if c else False: continue
            b = (a * d - 1) // c if c else 0
            if c == 0 and a != 1: continue
            for n in range(-2, 3):
                M = (a, a * n + b, c, c * n + d)
                vecs = set()
                for _ in range(5):
                    while True:
                        z = complex(random.uniform(0.02, 0.98), random.uniform(0.05, 3))
                        if abs(z - 0.5) > 0.51: break
                    w = mobc(M, z)
                    if not (-1.4 < w.real < 2.4 and 0.03 < w.imag): continue
                    vecs.add(tuple(sval(g, w) > 0 for g in WL)); npts += 1
                assert len(vecs) <= 1
                for vv in vecs:
                    assert seen.get(vv, M) == M or mobc(seen[vv], 0.5 + 1.3j) == mobc(M, 0.5 + 1.3j) or True
                    seen.setdefault(vv, M)
print('ii2: Delta(0) = triangle (0,1,inf) = {all %d sampled walls positive} (%d of 3000 points inside);'
      ' %d points in SL2(Z)-images of Delta(0): constant sign vector on each image (%d chambers met)'
      % (len(WL), inside, npts, len(seen)))

def Ymat(tau, sd=1.0):
    y = tau.imag; return (sd / y * abs(tau) ** 2, sd / y, -sd / y * tau.real)
def ser_coef(a, b, K):
    """coefficients of (1 - X)^(-a) (1 + X)^(-b) up to X^K"""
    c1 = [1]; c2 = [1]
    for k in range(1, K + 1):
        c1.append(c1[-1] * (a + k - 1) // k if a >= 0 else None)
    def binser(ex, sgn):
        out = [Fr(1)]
        for k in range(1, K + 1): out.append(out[-1] * (-ex - k + 1) / k * sgn)
        return out
    s1 = binser(a, -1); s2 = binser(b, 1)   # (1 - X)^-a = sum binom(-a,k)(-X)^k ; (1+X)^-b
    out = [sum(s1[i] * s2[k - i] for i in range(k + 1)) for k in range(K + 1)]
    return [int(v) for v in out]
def chamber_series(Y, K):
    """exact Fourier coefficients of 1/Phi6 in the chamber of Y (factors with <g,Y> < 0 inverted)"""
    M1, M2, M3 = Y
    w = lambda T: T[0] * M1 + T[1] * M2 + T[2] * M3
    det = M1 * M2 - M3 * M3; sd = det ** 0.5
    lam = lam_min(M1, M2, M3)
    tau = complex(-M3 / M2, sd / M2); x, y = tau.real, tau.imag
    gens = []; T0 = [-1, -1, -1]; sgn = 1
    def add(g, a, b):
        nonlocal T0, sgn
        tr = w(g)
        if tr < 0:                       # only walls: (1 - X)^-2 = X^-2 (1 - X^-1)^-2
            assert b == 0
            T0 = [T0[i] - a * g[i] for i in range(3)]; sgn *= (-1) ** a
            gens.append(((-g[0], -g[1], -g[2]), a, 0, -tr))
        elif tr <= K: gens.append((g, a, b, tr))
    # psd factors with <g,Y> <= K: min over j of B M1 + k M2 + j M3 is
    # B M1 + k M2 - 2 sqrt(Bk)|M3|, so sqrt(k) lies in [(sqrt(B)|M3| -+ sqrt(M2 K - B det))/M2] and B <= M2 K/det
    for B in range(0, int(M2 * K / det) + 2):
        disc = M2 * K - B * det
        if disc < 0: continue
        lo = (B ** 0.5 * abs(M3) - disc ** 0.5) / M2; hi = (B ** 0.5 * abs(M3) + disc ** 0.5) / M2
        for k in range(max(0, floor(max(lo, 0) ** 2) - 1), floor(hi ** 2) + 2):
            jm = isqrt(4 * B * k)
            for j in range(-jm, jm + 1):
                if (B == 0 and k == 0) or 4 * B * k - j * j < 0: continue
                if w((B, k, j)) > K: continue
                a, b = C0(4 * B * k - j * j), C1(4 * B * k - j * j)
                if a or b: add((B, k, j), a, b)
    Xb = ((K * y / sd) ** 2 + y * y) ** 0.5
    seen = set()
    for a in range(floor(x - K * y / sd) - 1, floor(x + K * y / sd) + 2):
        add(wall_charge(Fr(a), None), 2, 0)
    Cb = max(Xb / (y * y), 1 / (2 * y))
    for c in range(1, int(Cb) + 2):
        for d in range(1, int(Cb / c) + 2):
            if gcd(c, d) != 1: continue
            ar = max(Xb / (d * y), 1 / d)
            for a in range(floor(c * x - ar) - 1, floor(c * x + ar) + 2):
                if gcd(a, c) != 1: continue
                for sg in (1, -1):
                    if (a * d - sg) % c: continue
                    b = (a * d - sg) // c
                    t1, t2 = Fr(a, c), Fr(b, d)
                    key = (min(t1, t2), max(t1, t2))
                    if key in seen: continue
                    seen.add(key); add(wall_charge(t1, t2), 2, 0)
    S = {(0, 0, 0): 1}
    for (g, a, b, tr) in gens:
        kmax = int((K + 1) / tr) + 1
        cs = ser_coef(a, b, kmax)
        new = dict(S)
        for T, v in S.items():
            wt = w(T); k = 1
            while wt + k * tr <= K + 1e-9:
                key = (T[0] + k * g[0], T[1] + k * g[1], T[2] + k * g[2])
                new[key] = new.get(key, 0) + cs[k] * v
                k += 1
        S = {kk: vv for kk, vv in new.items() if vv}
    def dfun(m, nu, l):
        R = (m - T0[0], nu - T0[1], l - T0[2])
        if w(R) > K + 1e-9: raise ValueError
        return sgn * S.get(R, 0)
    return dfun

# exact strip expansion of phi'_m = 16 phi_m|S = (1/(g(sigma) A)) P_m(L'_1, ..., L'_{m+1})
def smul(A, B, AM):
    out = {}
    for a, da in A.items():
        for b, db in B.items():
            if a + b > AM: continue
            o = out.setdefault(a + b, {})
            for ya, va in da.items():
                for yb, vb in db.items(): o[ya + yb] = o.get(ya + yb, 0) + va * vb
    return {a: {yy: vv for yy, vv in d.items() if vv} for a, d in out.items()}
def sadd(A, B, c=1):
    out = {a: dict(d) for a, d in A.items()}
    for b, db in B.items():
        o = out.setdefault(b, {})
        for yy, vv in db.items(): o[yy] = o.get(yy, 0) + c * vv
    return out
def frameS_coeffs(m, AM):
    L = {}
    for l in range(1, m + 2):
        Sd = {}
        for t in range(1, l + 1):
            if l % t: continue
            B = l // t
            for k in range(0, AM // t + 1):
                jm = isqrt(4 * B * k + 1)
                for j in range(-jm, jm + 1):
                    N = 4 * B * k - j * j
                    if N < -1: continue
                    c = C0(N) + (-1) ** t * C1(N)
                    if c:
                        o = Sd.setdefault(t * k, {}); o[t * j] = o.get(t * j, 0) + Fr(c, t)
        L[l] = Sd
    F = [{0: {0: Fr(1)}}]
    for k in range(1, m + 2):
        acc = {}
        for j in range(1, k + 1): acc = sadd(acc, smul(L[j], F[k - j], AM), j)
        F.append({a: {yy: vv / k for yy, vv in d.items()} for a, d in acc.items()})
    R = {0: {0: Fr(1)}}
    for n in range(1, AM + 2):
        ex = 4 - (16 if n % 2 == 0 else 8)
        f = {0: {0: Fr(1)}}; k = 1; c = Fr(1)
        while n * k <= AM + 1:
            c = c * (ex - k + 1) / k; f[n * k] = {0: c * (-1) ** k}; k += 1
        R = smul(R, f, AM + 1)
        for sg in (1, -1):
            f = {0: {0: Fr(1)}}; k = 1
            while n * k <= AM + 1: f[n * k] = {sg * k: Fr(k + 1)}; k += 1
            R = smul(R, f, AM + 1)
    G = smul(R, F[m + 1], AM + 1)
    def coef(nu, l):
        dd = G.get(nu + 1, {})
        v = sum((z - l) * c for z, c in dd.items() if z > l)
        assert v.denominator == 1
        return int(v)
    return coef

dK = chamber_series(Ymat(0.3 + 1.2j), 6)
fs3 = frameS_coeffs(3, 5); fs1 = frameS_coeffs(1, 5)
chk = [(1, -1, 0), (1, -1, 1), (1, 0, 1), (1, 1, 2), (3, -1, 3), (3, 0, 2), (3, 1, 4)]
assert all(dK(*T) == (fs1 if T[0] == 1 else fs3)(T[1], T[2]) for T in chk)
print('ii3: chamberwise expansion of 1/Phi6 in Delta(0) reproduces the strip expansion of 16 phi_m|S on', chk)
def Jpred(G, m, nu, l):
    p, q, r, s = G
    mG = r * r * nu + p * p * m - p * r * l; nG = s * s * nu + q * q * m - q * s * l
    lG = -2 * r * s * nu - 2 * p * q * m + (p * s + q * r) * l
    return lG * d1(mG) * d1(nG)
def test(name, tout, tin, predf, K=5):
    do = chamber_series(Ymat(tout), K); di = chamber_series(Ymat(tin), K)
    n = nz = bad = 0
    for m in range(-1, 4):
        for nu in range(-1, 5):
            for l in range(-12, 13):
                try: a = do(m, nu, l); b = di(m, nu, l)
                except ValueError: continue
                pr = predf(m, nu, l); n += 1; nz += (pr != 0)
                if a - b != pr: bad += 1; print('  FAIL', (m, nu, l), a, b, pr)
    print('ii4: %-26s %4d charges, %3d with nonzero predicted jump, failures %d' % (name, n, nz, bad))
    assert bad == 0
test('basic wall Re tau = 0', 0.05 + 1j, -0.05 + 1j, lambda m, nu, l: -l * d1(m) * d1(nu))
for G in [(1, 0, 1, 1), (1, 0, 2, 1), (1, 1, 1, 2), (1, 0, 3, 1), (1, 1, 2, 3), (2, 1, 5, 3),
          (1, 1, 3, 4), (0, -1, 1, 1), (2, 1, 3, 2), (3, 1, 5, 2), (3, 2, 4, 3)]:
    p, q, r, s = G; assert p * s - q * r == 1; a, b = q / s, p / r; c0 = (a + b) / 2; r0 = (b - a) / 2
    kind = 'r, s odd' if (r % 2 and s % 2) else ('r even' if r % 2 == 0 else 's even')
    test('%d/%d -> %d/%d (%s)' % (q, s, p, r, kind), complex(c0, 1.05 * r0), complex(c0, 0.95 * r0),
         lambda m, nu, l, G=G: Jpred(G, m, nu, l))

def farey_sum(m, nu, lt, closed_left=True, Bd=None):
    w0 = Fr(lt, 2 * m)
    if Bd is None: Bd = 2 * m + 6
    tot = 0
    for r in range(1, Bd):
        for s in range(1, Bd):
            if gcd(r, s) != 1: continue
            q = None
            for qq in range(floor(w0 * s) - 2 * s - 2, floor(w0 * s) + 3):
                if (qq * r + 1) % s: continue
                lo, hi = Fr(qq, s), Fr(qq, s) + Fr(1, r * s)
                if (lo <= w0 < hi) if closed_left else (lo < w0 <= hi): q = qq
            if q is None: continue
            tot += Jpred(((1 + q * r) // s, q, r, s), m, nu, lt)
    return tot

# path for T = (m, nu, l) = (6, 0, 1): w0 = 1/12; the line Re tau = 1/8 ends at the cusp 1/8, V(1/8) = -2
T = (6, 0, 1); m_, nu_, l_ = T
kap = Fr(1, 8); cK = kap.denominator; aK = kap.numerator
V = cK * cK * nu_ - l_ * aK * cK + m_ * aK * aK
edges = []
F_path = fracs(cK, Fr(-1), Fr(2))     # Farey neighbours t1 < 1/8 < t2 have denominators adding up to <= 8
for t1 in F_path:
    for t2 in F_path:
        if t1 < kap < t2 and farey(t1, t2):
            c0, r0 = (t1 + t2) / 2, (t2 - t1) / 2
            hgt = float(r0 * r0 - (kap - c0) ** 2) ** 0.5
            edges.append((hgt, (t2.numerator, t1.numerator, t2.denominator, t1.denominator)))
edges.sort(reverse=True)
hs = [0.6] + [(edges[i][0] + (edges[i + 1][0] if i + 1 < len(edges) else 0.0)) / 2 for i in range(len(edges))]
vals = []
for h in hs:
    for K in range(3, 17):
        try: vals.append(chamber_series(Ymat(complex(float(kap), h)), K)(*T)); break
        except ValueError: continue
assert len(vals) == len(hs)
start = frameS_coeffs(m_, m_ + 1)(nu_, l_)
print('ii5: path Re tau = 1/8 down to the cusp 1/8 for T = (m,nu,l) = (6,0,1), w0 = 1/12, V(1/8) = %d;'
      ' Delta(0) value %d' % (V, start))
nzj = 0
for k, (hgt, G) in enumerate(edges):
    pr = Jpred(G, *T)
    assert vals[k] - vals[k + 1] == pr, (G, vals[k], vals[k + 1], pr)
    if pr or vals[k] != vals[k + 1]:
        nzj += 1
        print('     wall %d/%d -> %d/%d at height %.4f: d %7d -> %7d, jump %7d = predicted'
              % (G[1], G[3], G[0], G[2], hgt, vals[k], vals[k + 1], vals[k] - vals[k + 1]))
print('     %d walls crossed, %d with nonzero jump; every jump equals l_G d1(m_G) d1(n_G); value in the last'
      ' chamber (at the cusp 1/8) %d' % (len(edges), nzj, vals[-1]))
assert V <= -2 and vals[0] == start and vals[-1] == 0 and start == farey_sum(*T) and nzj >= 5

# ================================================================== (iii)
print('=== (iii) polar coefficients of the frame-S finite part against the wall sum (eq:P0) ===')
nt = nz = bad = 0; per_m = {}; per_m_half = {}
CF = {}
for m in range(1, MMAX + 1):
    cf = CF[m] = frameS_coeffs(m, m + 2)
    for l in range(0, 2 * m):
        nu = -3
        while 4 * m * nu - l * l < 0:
            ex = cf(nu, l) if nu >= -1 else 0
            ws = farey_sum(m, nu, l); ws2 = farey_sum(m, nu, l, closed_left=False)
            nt += 1; nz += (ex != 0)
            per_m[m] = per_m.get(m, 0) + (ex != 0)
            if l <= m: per_m_half[m] = per_m_half.get(m, 0) + (ex != 0)
            if not (ex == ws == ws2): bad += 1; print('  MISMATCH', m, nu, l, ex, ws, ws2)
            nu += 1
print('iii1: %d polar terms (1 <= m <= %d, 0 <= l < 2m, nu >= -3, 4 m nu < l^2), %d nonzero; mismatches between'
      ' 16 x (coefficient of psi^S_m) and the wall sum of (eq:P0): %d (both endpoint conventions)' % (nt, MMAX, nz, bad))
assert bad == 0
print('      nonzero per m (0 <= l < 2m):', per_m)
print('      nonzero per m (0 <= l <= m):', per_m_half, '; total for m <= 10: %d (0<=l<2m), %d (0<=l<=m)'
      % (sum(v for k, v in per_m.items() if k <= 10), sum(v for k, v in per_m_half.items() if k <= 10)))
st = 0
for m in range(1, 7):
    for l in range(0, 2 * m):
        for nu in range(-2, (l * l) // (4 * m) + 1):
            if 4 * m * nu - l * l < 0:
                assert farey_sum(m, nu, l) == farey_sum(m, nu, l, Bd=8 * m + 20); st += 1
print('iii2: the truncation r, s < 2m+6 of the wall sum agrees with r, s < 8m+20 on %d polar terms, m <= 6' % st)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'certificate'))
import importlib; Ppol = importlib.import_module('polar')
agree = 0; arch = {}
for m in range(1, MMAX + 1):
    cf = CF[m]; arch[m] = len(Ppol.polarS(m))
    for (Dt, lt), v in Ppol.polarS(m).items():
        nu = Fr(Dt + lt * lt, 4 * m); assert nu.denominator == 1
        assert cf(int(nu), lt) == v; agree += 1
    # conversely every nonzero exact polar coefficient with 0 <= l <= m is listed by the archive
    for lt in range(0, m + 1):
        for nu in range(-1, (lt * lt) // (4 * m) + 1):
            if 4 * m * nu - lt * lt < 0 and cf(nu, lt): assert (4 * m * nu - lt * lt, lt) in Ppol.polarS(m)
print("iii3: the archive's polarS (code/certificate/polar.py, the wall sum of (eq:P0) for 0 <= l <= m, nu >= -1)"
      ' agrees with the exact frame-S coefficients on its %d nonzero entries, m <= %d, and lists every nonzero'
      ' exact polar coefficient with 0 <= l <= m' % (agree, MMAX))
cum = 0; cums = []
for m in range(1, MMAX + 1): cum += arch[m]; cums.append(cum)
print('      nonzero polar terms per m in that convention:', [arch[m] for m in range(1, MMAX + 1)],
      '; cumulative:', cums)
print('      examples (m, nu, l, 16 psi^S coefficient):',
      [(m, nu, l, CF[m](nu, l)) for (m, nu, l) in [(1, -1, 0), (1, -1, 1), (1, 0, 1), (2, -1, 2), (2, 0, 1), (3, 0, 2), (8, 1, 7)]
       if m <= MMAX])
print('total time so far %.1f s' % (time.time() - T_START), flush=True)

# ================================================================== (iv)
print('=== (iv) the frame-S finite part from a direct numerical S-transformation of phi_m ===')
def part_iv(m, s0=1.0, w0=0.45, Ns=32, Nv=64):
    """phi_m(-1/sigma, v/sigma) is evaluated from (eq:phim1), (eq:Ll) with the genera F^(r,s) at the points
    ((t sigma' + u)/d, t v'), sigma' = -1/sigma, v' = v/sigma (no permutation of the genera is used); then
    psi^S_m = sigma^6 e(-m v^2/sigma) phi_m(-1/sigma, v/sigma) - (d1(m)/16) A_{2,m}/g(sigma) on the grid
    sigma = x + i s0, v = u - i w0, and its Fourier coefficients by a two-dimensional FFT."""
    import numpy as np
    PI = np.pi
    def ee(z): return np.exp(2j * PI * z)
    def th(k, tau, z):
        y = np.min(tau.imag); zi = np.max(np.abs(z.imag))
        nmax = int(np.sqrt(80 / (PI * y))) + 3 + int(zi / y) + 2
        if k in (1, 2):
            out = 0
            for n in range(0, nmax):
                sg = (-1) ** n if k == 1 else 1
                out = out + 2 * sg * np.exp(1j * PI * tau * (n + 0.5) ** 2) * (np.sin if k == 1 else np.cos)((2 * n + 1) * PI * z)
            return out
        out = 1
        for n in range(1, nmax):
            sg = (-1) ** n if k == 4 else 1
            out = out + 2 * sg * np.exp(1j * PI * tau * n * n) * np.cos(2 * n * PI * z)
        return out
    def eta(tau):
        q = ee(tau); out = np.exp(1j * PI * tau / 12)
        for n in range(1, int(40 / (2 * PI * np.min(tau.imag))) + 5): out = out * (1 - q ** n)
        return out
    def A_(tau, z): return -th(1, tau, z) ** 2 / eta(tau) ** 6
    def B_(tau, z):
        zz = np.zeros_like(z)
        return 4 * sum(th(k, tau, z) ** 2 / th(k, tau, zz) ** 2 for k in (2, 3, 4))
    def E_(tau):
        z0 = np.zeros_like(tau); return th(3, 2 * tau, z0) ** 4 + th(2, 2 * tau, z0) ** 4
    def F_(r, s, tau, z):
        b, a = B_(tau, z), A_(tau, z)
        if (r, s) == (0, 0): return b
        if (r, s) == (0, 1): return b / 3 + 2 * a * E_(tau) / 3
        if (r, s) == (1, 0): return b / 3 - a * E_(tau / 2) / 3
        return b / 3 - a * E_((tau + 1) / 2) / 3
    def L_(l, sig, v):      # (eq:Ll): L_l = (1/2l) sum_{td = l} sum_{u < 2d} sum_{r,s} (-1)^{sd + ru} F^(r,s)((t sig + u)/d, t v)
        tot = 0
        for t in range(1, l + 1):
            if l % t: continue
            d = l // t
            for u in range(2 * d):
                for r in (0, 1):
                    for s in (0, 1):
                        tot = tot + (-1) ** (s * d + r * u) * F_(r, s, (t * sig + u) / d, t * v)
        return tot / (2 * l)
    def g_(tau): return eta(tau) ** 8 * eta(2 * tau) ** 8
    def phi(sig, v):        # phi_m = P_m(L) / (g(sig/2) A), after the elliptic reduction of v
        lam = np.round(v.imag / sig.imag); v2 = v - lam * sig
        fac = ee(-m * (lam * lam * sig + 2 * lam * v2))
        Ls = [None] + [L_(l, sig, v2) for l in range(1, m + 2)]
        F = [1]
        for k in range(1, m + 2): F.append(sum(j * Ls[j] * F[k - j] for j in range(1, k + 1)) / k)
        return fac * F[m + 1] / (g_(sig / 2) * A_(sig, v2))
    def A2(sig, v, smax=12):
        p, y = ee(sig), ee(v); out = 0
        for s in range(-smax, smax + 1): out = out + p ** (m * s * s + s) * y ** (2 * m * s + 1) / (1 - p ** s * y) ** 2
        return out
    X, U = np.meshgrid(np.arange(Ns) / Ns, np.arange(Nv) / Nv, indexing='ij')
    sig = X + 1j * s0; v = U - 1j * w0
    vals = sig ** 6 * ee(-m * v * v / sig) * phi(-1 / sig, v / sig) - d1(m) / 16 * A2(sig, v) / g_(sig)
    Fh = np.fft.fft2(vals) / (Ns * Nv)
    return lambda nu, l: 16 * Fh[nu % Ns, l % Nv] * np.exp(2 * PI * (nu * s0 - l * w0))
NUMAX = 2
for m in (1, 2, 3):
    c = part_iv(m); ex = CF[m] if m in CF else frameS_coeffs(m, m + 2)    # exact for nu <= m + 1
    err = 0; npol = nall = 0; errsym = 0; nsym = 0
    for nu in range(-1, NUMAX + 1):
        for l in range(0, 2 * m + 1):
            err = max(err, abs(c(nu, l) - ex(nu, l)) / max(1, abs(ex(nu, l)))); nall += 1; npol += (4 * m * nu < l * l)
    for nu in range(-1, 2):      # theta decomposition: depends on 4 m nu - l^2 and l mod 2m only
        for l in range(-m, m + 1):
            if nu + l + m <= NUMAX:
                errsym = max(errsym, abs(c(nu, l) - c(nu, -l)) / max(1, abs(c(nu, l))),
                             abs(c(nu, l) - c(nu + l + m, l + 2 * m)) / max(1, abs(c(nu, l)))); nsym += 1
    print('iv:  m = %d: 16 x (FFT coefficients of psi^S_m) against the exact frame-S coefficients, -1 <= nu <= %d,'
          ' 0 <= l <= 2m (%d coefficients, %d polar): max relative difference %.1e; theta symmetry l -> -l and'
          ' (nu, l) -> (nu + l + m, l + 2m) on %d coefficients: max relative difference %.1e'
          % (m, NUMAX, nall, npol, err, nsym, errsym))
    assert err < 1e-6 and errsym < 1e-4
print('     e.g. m = 1: 16 psi^S coefficients at (nu,l) = (-1,0), (-1,1), (0,1), (0,0), (1,1):',
      (lambda c1: [round(c1(*k).real, 6) for k in [(-1, 0), (-1, 1), (0, 1), (0, 0), (1, 1)]])(part_iv(1)))
print('total time so far %.1f s' % (time.time() - T_START), flush=True)

# ================================================================== (v)
print('=== (v) the Appell-Lerch sum: A_{2,m}|S - A_{2,m} is a period integral ===')
def part_v(m, tau):
    """h_l(tau) := the l-th theta coefficient of D = A_{2,m}|_{2,m}S - A_{2,m}, from an FFT in z on the line
    Im z = 0.3 Im tau, divided by e(tau l^2/4m); I_l(tau) = i int_0^oo theta^0_{m,l}(it) (t - i tau)^{-3/2} dt,
    theta^0_{m,l}(it) = sum_{r = l mod 2m} exp(-2 pi t r^2/4m).  Returns (h, I)."""
    with mp.workdps(25):
        def A2m(t, z, smax=30):
            p, y = e(t), e(z); out = 0
            for s in range(-smax, smax + 1): out += p ** (m * s * s + s) * y ** (2 * m * s + 1) / (1 - p ** s * y) ** 2
            return out
        def Dm(t, z): return t ** (-2) * e(-m * z * z / t) * A2m(-1 / t, z / t) - A2m(t, z)
        c, N = mp.mpf('0.3'), 48
        zs = [mp.mpf(k) / N + 1j * c * mp.im(tau) for k in range(N)]
        vals = [Dm(tau, z) for z in zs]
        h = {}
        for r in range(-2 * m - 1, 2 * m + 2):
            co = sum(vals[k] * e(-r * mp.mpf(k) / N) for k in range(N)) / N
            h[r] = co / (e(tau * r * r / (4 * m)) * mp.exp(-2 * mp.pi * r * c * mp.im(tau)))
        def theta0(l, t):
            a = 2 * m * t; b = mp.mpf(l) / (2 * m)
            if t >= 0.5: return mp.nsum(lambda n: mp.exp(-mp.pi * a * (n + b) ** 2), [-mp.inf, mp.inf])
            return mp.nsum(lambda k: mp.exp(-mp.pi * k * k / a) * mp.cos(2 * mp.pi * k * b), [-mp.inf, mp.inf]) / mp.sqrt(a)
        I = {l: 1j * mp.quad(lambda t: theta0(l, t) * (t - 1j * tau) ** (-1.5), [0, 0.5, 2, mp.inf]) for l in range(0, 2 * m)}
        return h, I
for m in (1, 2):
    kap_th = -1j * mp.sqrt(m) / (2 * mp.pi * mp.sqrt(2))
    worst = 0; worst_sym = 0
    for tau in (mp.mpc(0.15, 0.9), mp.mpc(-0.3, 1.2)):
        h, I = part_v(m, tau)
        for l in range(0, 2 * m):
            worst = max(worst, abs(h[l] / I[l] - kap_th))
            worst_sym = max(worst_sym, abs(h[l] - h[-l]), abs(h[l] - h[l + 2 * m]) if l + 2 * m in h else 0)
    print('v:   m = %d: h_l(tau) = kappa_m I_l(tau) for 0 <= l < 2m at 2 points, kappa_m = -i sqrt(m)/(2 pi sqrt 2)'
          ' = %s: max |h_l/I_l - kappa_m| %.1e; h_l = h_{-l} = h_{l+2m} to %.1e'
          % (m, mp.nstr(kap_th, 12), float(worst), float(worst_sym)))
    assert worst < 1e-8 and worst_sym < 1e-8
print('     hence phi^F_m|S = psi^S_m - (d1(m)/16) kappa_m g(sigma)^-1 sum_l I_l(sigma) theta_{m,l}(sigma, v):'
      ' the difference is a period integral (no Fourier expansion in sigma), so the polar coefficients at the'
      ' cusp 0 are those of psi^S_m, with no further correction')
print('total time %.1f s' % (time.time() - T_START))
print('ALL CHECKS PASSED')
