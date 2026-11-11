"""Exact Fourier coefficients of the Fourier-Jacobi coefficients of 1/Phi~6 (frame I, cusp infinity)
and of phi'_m/16 = P_m(L'_1..L'_{m+1})/(16 g(sigma) phi_{-2,1}) (frame S), with L'_l the series of polar0
Lemma lem:LlS; the latter are also the coefficients of -1/(16 Phi6^DJS)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from fractions import Fraction as Fr
from math import isqrt
import flint
import chl_coefficients as C

def frameI(Lq, Ku):
    coef, _ = C.chamber_coefficients(2, Lq, Ku)
    return coef          # coef(m, n, l): coefficient of e(m rho) x^n y^l, x = e(sigma/2), |y| > 1 strip

def exps_S(Dmax):
    """C0(D) = sum_s c^(0,s)(D), C1(D) = sum_s c^(1,s)(D) from F^(0,0)+F^(0,1) = 4/3 B + 2/3 A E,
    F^(1,0)+F^(1,1) = 2/3 (B - A E), E = 2E2(2tau) - E2(tau)."""
    A, Bm = C.disc_tables(Dmax)          # A[D] = coeff of B/4, Bm[D] = coeff of A (phi_{-2,1})
    def cB(D): return 4 * A.get(D, 0) if D >= -1 else 0
    def cA(D): return Bm.get(D, 0) if D >= -1 else 0
    T = Dmax // 4 + 2
    def s1odd(n): return sum(d for d in range(1, n + 1) if n % d == 0 and d % 2)
    e = [1] + [24 * s1odd(t) for t in range(1, T + 1)]
    C0, C1 = {}, {}
    for D in range(-1, Dmax + 1):
        EA = sum(e[t] * cA(D - 4 * t) for t in range(0, (D + 1) // 4 + 1))
        c0 = Fr(4, 3) * cB(D) + Fr(2, 3) * EA
        c1 = Fr(2, 3) * (cB(D) - EA)
        assert c0.denominator == 1 and c1.denominator == 1
        C0[D] = int(c0); C1[D] = int(c1)
    return C0, C1

def frameS(Lq, Ku):
    """coefficient of e(m rho) p^n y^l (p = e(sigma)) of -1/(16 Phi6^DJS), expanded with |y| > 1,
    Phi6 = -e(sigma+rho+v) prod_{(k,l,b)>0} (1 - (-1)^r e(k sigma + l rho + b v))^{c^(r,s)(4kl-b^2)}."""
    C0, C1 = exps_S(4 * Lq * Ku + 8)
    R = Lq + Ku; W = 4 * R + 3
    def enc(d):
        L = [0] * ((Ku + 1) * W)
        for (b, r), v in d.items(): L[b * W + r + R] = flint.fmpq(v.numerator, v.denominator)
        return flint.fmpq_poly(L)
    def dec(p):
        out = {}
        for i, v in enumerate(p.coeffs()):
            if v:
                b, rr = divmod(i, W); r = rr - 2 * R
                if b <= Ku:
                    out[(b, r)] = Fr(int(v.p), int(v.q))
        return out
    Lg = [dict() for _ in range(Lq + 1)]
    for l in range(Lq + 1):
        for u in range(Ku + 1):
            if l == u == 0: continue
            jm = isqrt(4 * u * l + 1) + 1
            for r in range(-jm, jm + 1):
                D = 4 * u * l - r * r
                if D < -1: continue
                c0, c1 = C0.get(D, 0), C1.get(D, 0)
                if not (c0 or c1): continue
                t = 1
                while t * l <= Lq and t * u <= Ku:
                    key = (t * u, t * r)
                    Lg[t * l][key] = Lg[t * l].get(key, 0) + Fr(c0 + (-1) ** t * c1, t)
                    t += 1
    assert C0[-1] == 2 and C1[-1] == 0
    L0 = {}
    for (b, r), v in Lg[0].items(): L0.setdefault(b, {})[r] = v
    F0 = {0: {0: Fr(1)}}
    for b in range(1, Ku + 1):
        acc = {}
        for j in range(1, b + 1):
            for r1, v1 in L0.get(j, {}).items():
                for r2, v2 in F0[b - j].items():
                    acc[r1 + r2] = acc.get(r1 + r2, 0) + j * v1 * v2
        F0[b] = {r: v / b for r, v in acc.items() if v}
    F = [{(b, r): v for b, d in F0.items() for r, v in d.items()}]
    Lp = [enc(d) for d in Lg]; Fp = [enc(F[0])]
    for k in range(1, Lq + 1):
        acc = flint.fmpq_poly(0)
        for j in range(1, k + 1):
            if Lg[j]: acc += j * Lp[j] * Fp[k - j]
        d = {key: v / k for key, v in dec(acc).items()}
        F.append(d); Fp.append(enc(d))
    rows = {}
    for a, d in enumerate(F):
        for (b, z), v in d.items():
            if v: rows.setdefault((a, b), {})[z] = v
    def coef(m, n, l):
        assert -1 <= m < Lq and -1 <= n < Ku
        v = sum((z - l) * c for z, c in rows.get((m + 1, n + 1), {}).items() if z > l)
        return Fr(v, 16)
    return coef

if __name__ == '__main__':
    cs = frameS(3, 12)
    # psi^S_{-1} = 1/(16 g(sigma) A): A = y - 2 + y^-1 + O(p), 1/g = p^-1 + 8 + ...
    for l in range(-3, 4):
        print(l, [cs(-1, n, l) for n in range(-1, 5)])
    print('m=1', [[cs(1, n, l) for n in range(-1, 6)] for l in range(0, 3)])
