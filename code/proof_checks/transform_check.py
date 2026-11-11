#!/usr/bin/env python3
"""Numerical checks for transform.tex (the transformation law (T2b) and its ingredients).

Checks, for m = 1, 2, 3 (mpmath, 50 digits):
 (a) Lemma tr:lem:theta: theta_{m,l'}(M tau, z/J) = J^(1/2) e(m g z^2/J) sum_l W_{g,d}(l',l) theta_{m,l}(tau,z);
 (b) Lemma tr:lem:Afourier: the Fourier expansion of A_{2,m} in the strip -Im tau < Im z < 0;
 (c) Proposition tr:prop:Ahat: A^ = A_{2,m} - kappa_m sum_l E_l theta_{m,l} satisfies A^|_{2,m} S = A^, A^|T = A^;
     and the period identity of its proof, (A|S - A) = sum_l c_l theta_{m,l} with
     c_l = (sqrt(2m)/4pi) e^{i pi/4} int_0^{i inf} theta_l(w)(w+tau)^(-3/2) dw;
 (d) Proposition tr:prop:jacobi(iii), i.e. polar0 Proposition prop:PhiS(iii): psi_m|_{-6,m} S = phi'_m/16, with
     the coefficients of psi_m = phi_{-1} P_m(L) computed exactly from the product formula of Phi~6 (cusp infinity),
     and those of phi'_m = P_m(L'_1..L'_{m+1})/(g(sigma) phi_{-2,1}) computed exactly from the series
     L'_l = sum_t (1/t) sum (C0 + (-1)^t C1)(4(l/t)k - j^2) p^{tk} y^{tj} of polar0 Lemma lem:LlS (tc_series.frameS);
     both sides evaluated from finite parts (theta series) plus G_F A_{2,m};
 (e) (T2a) and (T2b) for several (gamma, delta) in all three frames.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from math import gcd
from tc_series import frameI, frameS

mp.mp.dps = 50
I = mp.mpc(0, 1)
def e(x): return mp.exp(2 * mp.pi * I * x)

KU = 200
cI = frameI(5, KU)
cS = frameS(5, KU)

def dvals(K):
    """coefficients d1(k), k >= -1, of 1/(eta^8(tau) eta^8(2 tau))"""
    N = K + 3
    P = [0] * (N + 1); P[0] = 1
    for n in range(1, N + 1):
        for e_ in ([8] if n % 2 else [8]) :
            pass
    # prod (1-q^n)^8 (1-q^{2n})^8
    P = [0] * (N + 1); P[0] = 1
    def mulfac(P, step, pw):
        for _ in range(pw):
            for k in range(N, step - 1, -1): P[k] -= P[k - step]
    for n in range(1, N + 1):
        mulfac(P, n, 8)
        if 2 * n <= N: mulfac(P, 2 * n, 8)
    inv = [0] * (N + 1); inv[0] = 1
    for n in range(1, N + 1): inv[n] = -sum(P[k] * inv[n - k] for k in range(1, n + 1))
    return inv                       # 1/g = q^-1 sum inv[j] q^j
D1 = dvals(60)
def fhat(m): return D1[m + 1]

def g_(tau):                          # eta(tau)^8 eta(2 tau)^8
    q = e(tau); s = q
    for n in range(1, 400):
        t = (1 - q ** n) ** 8 * (1 - q ** (2 * n)) ** 8
        s *= t
        if abs(q ** n) < mp.mpf(10) ** (-60): break
    return s

def theta_null(m, l, tau):
    s = 0; r0 = l % (2 * m)
    R = int(mp.sqrt(4 * m * 140 / (2 * mp.pi * mp.im(tau)))) + 3
    for k in range(-R, R + 1):
        r = r0 + 2 * m * k; s += e(r * r * tau / (4 * m))
    return s

def theta_mz(m, l, tau, z):
    s = 0; r0 = l % (2 * m)
    R = int(mp.sqrt(4 * m * 140 / (2 * mp.pi * mp.im(tau)))) + int(abs(mp.im(z)) / mp.im(tau) * 2) + 6
    for k in range(-R, R + 1):
        r = r0 + 2 * m * k; s += e(r * r * tau / (4 * m) + r * z)
    return s

def W(m, g, d, lt, l):
    al = pow(d % g, -1, g) if g > 1 else 0
    s = 0
    for T in range(g):
        x = lt - 2 * m * T
        s += e(mp.mpf(al * x * x - 2 * l * x + d * l * l) / (4 * m * g))
    return s / mp.sqrt(2 * m * g * I)

def E_(m, l, tau):
    """E_l(tau) = int_{-conj tau}^{i inf} theta_l(w)(w+tau)^(-3/2) dw, closed form"""
    t2 = mp.im(tau); s = 0; r0 = l % (2 * m)
    R = int(mp.sqrt(4 * m * 140 / (2 * mp.pi * t2))) + 3
    for k in range(-R, R + 1):
        r = r0 + 2 * m * k
        if r == 0: s += 2 / mp.sqrt(2 * t2)
        else:
            a = mp.pi * r * r / (2 * m)
            s += e(-r * r * tau / (4 * m)) * mp.sqrt(a) * mp.gammainc(-0.5, 2 * a * t2)
    return mp.exp(-I * mp.pi / 4) * s

def kappa(m): return -mp.exp(I * mp.pi / 4) * mp.sqrt(m) / (2 * mp.sqrt(2) * mp.pi)

def A2(m, tau, z):
    q = e(tau); y = e(z); s = 0
    for k in range(-60, 61):
        s += q ** (m * k * k + k) * y ** (2 * m * k + 1) / (1 - q ** k * y) ** 2
    return s

def Ahat(m, tau, z):
    return A2(m, tau, z) - kappa(m) * sum(E_(m, l, tau) * theta_mz(m, l, tau, z) for l in range(2 * m))

def rel(a, b): return abs(a - b) / max(abs(a), abs(b), mp.mpf(10) ** -40)

out = []
def rep(s):
    print(s, flush=True); out.append(s)

worst = {}
def note(key, v):
    worst[key] = max(worst.get(key, 0), v)

# ------------------------------------------------------------------ (a) theta transformation
pts = [(mp.mpc('0.17', '0.83'), mp.mpc('0.31', '-0.22')), (mp.mpc('-0.4', '1.3'), mp.mpc('-0.12', '0.4'))]
for m in (1, 2, 3):
    for g in range(1, 5):
        for md in range(0, 2 * g):
            d = -md
            if gcd(g, d) != 1: continue
            al = pow(d % g, -1, g) if g > 1 else 0; be = (al * d - 1) // g
            for tau, z in pts:
                J = g * tau + d; Mt = (al * tau + be) / J
                for lp in range(2 * m):
                    lhs = theta_mz(m, lp, Mt, z / J)
                    rhs = mp.sqrt(J) * e(m * g * z * z / J) * sum(W(m, g, d, lp, l) * theta_mz(m, l, tau, z) for l in range(2 * m))
                    note('(a) theta lemma', rel(lhs, rhs))
rep('(a) theta transformation, m<=3, gamma<=4, all delta, l: max rel. error %s' % mp.nstr(worst['(a) theta lemma'], 3))

# ------------------------------------------------------------------ (b) Fourier expansion of A in the strip
def a_coef(m, r, tau):
    s = 0
    if r > 0:
        n = r - 2 * m
        while n >= 1: s += n * e(-n * n * tau / (4 * m)); n -= 2 * m
    elif r < 0:
        n = r
        while n <= -1: s += -n * e(-n * n * tau / (4 * m)); n += 2 * m
    return s * e(r * r * tau / (4 * m))
for m in (1, 2, 3):
    tau = mp.mpc('0.21', '0.9'); z = mp.mpc('0.37', '-0.41')
    lhs = A2(m, tau, z); rhs = sum(a_coef(m, r, tau) * e(r * z) for r in range(-160, 161))
    note('(b) A strip', rel(lhs, rhs))
rep('(b) Fourier expansion of A_{2,m} in the strip, m<=3: max rel. error %s' % mp.nstr(worst['(b) A strip'], 3))

# ------------------------------------------------------------------ (c) completion
def theta_it(m, l, t):     # theta_{m,l}(i t), t > 0; for t < 1 through theta_{m,l}(-1/tau') with tau' = i/t
    if t >= 1: return theta_null(m, l, I * t)
    return sum(e(-mp.mpf(l * lp) / (2 * m)) * theta_null(m, lp, I / t) for lp in range(2 * m)) / mp.sqrt(2 * m * t)
def Ptilde(m, l, tau):     # int_0^{i inf} theta_l(w)(w+tau)^(-3/2) dw
    f = lambda t: theta_it(m, l, t) * (I * t + tau) ** mp.mpf(-1.5) * I
    return mp.quad(f, [0, mp.mpf('0.05'), 1, 10, mp.inf])
for m in (1, 2, 3):
    for tau, z in [(mp.mpc('0.23', '0.97'), mp.mpc('0.31', '-0.27')), (mp.mpc('-0.35', '1.21'), mp.mpc('0.2', '0.33'))]:
        lhsS = tau ** -2 * e(-m * z * z / tau) * Ahat(m, -1 / tau, z / tau)
        note('(c) Ahat|S', rel(lhsS, Ahat(m, tau, z)))
        note('(c) Ahat|T', rel(Ahat(m, tau + 1, z), Ahat(m, tau, z)))
    # period identity at one point: (A|S - A)(tau, z) = sum_l c_l theta_{m,l}
    tau, z = mp.mpc('0.11', '0.93'), mp.mpc('0.27', '-0.31')
    lhs = tau ** -2 * e(-m * z * z / tau) * A2(m, -1 / tau, z / tau) - A2(m, tau, z)
    cst = mp.sqrt(2 * m) / (4 * mp.pi) * mp.exp(I * mp.pi / 4)
    rhs = sum(cst * Ptilde(m, l, tau) * theta_mz(m, l, tau, z) for l in range(2 * m))
    note('(c) period identity', rel(lhs, rhs))
    # closed form of E_l against quadrature
    tau = mp.mpc('0.3', '0.7')
    for l in range(2 * m):
        q = mp.quad(lambda t: theta_null(m, l, -mp.conj(tau) + I * t) * (-mp.conj(tau) + I * t + tau) ** mp.mpf(-1.5) * I, [0, 1, 10, mp.inf])
        note('(c) E closed form', rel(q, E_(m, l, tau)))
for k in ('(c) Ahat|S', '(c) Ahat|T', '(c) period identity', '(c) E closed form'):
    rep('%s, m<=3: max rel. error %s' % (k, mp.nstr(worst[k], 3)))

# ------------------------------------------------------------------ frame data
def hI(m, l, tau):
    return sum(cI(m, n, l) * e(mp.mpf(2 * m * n - l * l) * tau / (4 * m)) for n in range(-1, KU))
def hS(m, l, tau):
    return sum(mp.mpf(cS(m, n, l).numerator) / cS(m, n, l).denominator * e(mp.mpf(4 * m * n - l * l) * tau / (4 * m)) for n in range(-1, KU))
def hT(m, l, tau): return e(mp.mpf(l * l) / (4 * m)) * hI(m, l, tau + 1)
def GI(m, tau): return fhat(m) / g_(tau / 2)
def GT(m, tau): return GI(m, tau + 1)
def GS(m, tau): return fhat(m) / (16 * g_(tau))
H = {'I': hI, 'T': hT, 'S': hS}; G = {'I': GI, 'T': GT, 'S': GS}
def hhat(F, m, l, tau): return H[F](m, l, tau) + kappa(m) * G[F](m, tau) * E_(m, l, tau)

def frame(N):     # N = (a, b, c, d); F with N R_F^{-1} in Gamma^0(2)
    a, b, c, d = N
    if b % 2 == 0: return 'I'
    # N R_T^{-1} = (a, b - a, c, d - c)
    if (b - a) % 2 == 0: return 'T'
    return 'S'

# ------------------------------------------------------------------ (d) frame identification, uncompleted
for m in (1, 2, 3):
    sig = mp.mpc('0.02', '1.0'); v = mp.mpc('0.35', '-0.45')     # -Im sig < Im v < 0, and (sig', v') = (-1/sig, v/sig) in the strip
    s2, v2 = -1 / sig, v / sig
    assert -mp.im(s2) < mp.im(v2) < 0
    psiI = sum(hI(m, l, s2) * theta_mz(m, l, s2, v2) for l in range(2 * m)) + GI(m, s2) * A2(m, s2, v2)
    lhs = sig ** 6 * e(-m * v * v / sig) * psiI
    rhs = sum(hS(m, l, sig) * theta_mz(m, l, sig, v) for l in range(2 * m)) + GS(m, sig) * A2(m, sig, v)
    note('(d) psi_I|S = psi_S', rel(lhs, rhs))
rep('(d) psi_m|_{-6,m}S = phi_m\'/16 (series of L\'_l), m<=3: rel. error %s' % mp.nstr(worst['(d) psi_I|S = psi_S'], 3))

# ------------------------------------------------------------------ (e) (T2a), (T2b)
# tau chosen with |j| = 1, so that Im tau = Im N tau = (1 - 0.01)^(1/2)/gamma (the q-series of both frames converge
# at the truncation KU); gamma = 1, 2 meet all three frames, gamma = 3 the frames T and S.
cases = []
for g, d in [(1, 0), (1, -1), (2, -1), (2, -3), (3, -1), (3, -2), (3, -4)]:
    al = pow(d % g, -1, g) if g > 1 else 0
    cases.append((g, d, mp.mpf(al) / g + mp.mpc('0.1', mp.sqrt(mp.mpf('0.99'))) / g))
for m in (1, 2, 3):
    for g, d, tau in cases:
        al = pow(d % g, -1, g) if g > 1 else 0; be = (al * d - 1) // g
        N = (d, -be, -g, al); F = frame(N)
        j = -g * tau + al; Nt = (d * tau - be) / j
        e2a = rel(GI(m, Nt), j ** -8 * G[F](m, tau))
        note('(e) T2a', e2a)
        errs = []
        for l in range(2 * m):
            lhs = hhat('I', m, l, Nt)
            rhs = j ** mp.mpf(-6.5) * sum(W(m, g, d, lp, l) * hhat(F, m, lp, tau) for lp in range(2 * m))
            errs.append(rel(lhs, rhs))
        note('(e) T2b', max(errs))
        rep('(e) m=%d gamma=%d delta=%d frame %s, Im tau=%s, Im N tau=%s: (T2a) %s, (T2b) max over l %s'
            % (m, g, d, F, mp.nstr(mp.im(tau), 3), mp.nstr(mp.im(Nt), 3), mp.nstr(e2a, 3), mp.nstr(max(errs), 3)))
# negative control: drop the shadow term
m, g, d, tau = 2, 1, 0, mp.mpc('0.1', mp.sqrt(mp.mpf('0.99')))
j = -tau; Nt = -1 / tau
lhs = hI(m, 0, Nt); rhs = j ** mp.mpf(-6.5) * sum(W(m, 1, 0, lp, 0) * hS(m, lp, tau) for lp in range(2 * m))
rep('control: (T2b) without the completion terms, m=2, gamma=1: rel. error %s (must be large)' % mp.nstr(rel(lhs, rhs), 3))
rep('SUMMARY ' + '; '.join('%s %s' % (k, mp.nstr(v, 3)) for k, v in worst.items()))
