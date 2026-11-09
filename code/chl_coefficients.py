#!/usr/bin/env python3
"""Positivity of the single-centred index of the Z_N CHL models, N = 2, 3, 5, 7, for every
charge of bounded discriminant.

Charges: m = N Q^2/2, n = P^2/2, l = Q.P, Delta = 4mn - N l^2 = N (Q^2 P^2 - (Q.P)^2).
The duality group Gamma_1(N) acts by (Q, P) -> (aQ + bP, cQ + dP); on the charge form
f = m s^2 + N l s t + N n t^2, of discriminant -N Delta, this is f -> f o g, (f o g)(s, t) = f((s, t) g).

Steps, all in exact arithmetic:
 B1  the Gamma_1(N)-orbits of charges with given Delta, through Gauss reduction:
     an orbit is a reduced form F of discriminant -N Delta with a class of bottom rows mod N
     modulo the automorphs of F;
 B2  in every orbit a charge of the set A of Bhand-Sen-Singh, (5.46)-(5.47), with 0 <= l < 2m;
     the test is (5.47) with j = 2s + 1:  N (j m + l m1)^2 + m1^2 Delta >= N m^2 for every
     m1 in N Z, m1 != 0, and odd j with m1 | (j^2 - 1)/4;
 B3  the coefficient c(m, n, l) of q^m x^n y^l of 1/Phi~ in chamber R, from the product formula
     (exp-log recursion); the exponents are symmetric under q <-> x on the whole range used, so
     c(m, n, l) = c(n, m, l) there, and each charge is read in the cheaper orientation;
 B4  d* = (-1)^(l+1) c(m, n, l) is positive in every orbit, and agrees for
     all representatives of an orbit found in the box; at m = 1 it reproduces the first known values;
 B5  on the whole range, every charge of least m in its orbit lies in A.

Usage: python3 chl_coefficients.py [N:Dmax ...] [--out FILE] [--values FILE]
Default: 2:200 3:200 5:20 7:20.
"""
import sys, json, time, hashlib
from fractions import Fraction as Fr
from math import isqrt, gcd
import flint

# ---------------------------------------------------------------- B3: product exponents
def _conv(a, b, X):
    co = [int(v) for v in (flint.fmpz_poly(a[:X + 1]) * flint.fmpz_poly(b[:X + 1])).coeffs()]
    return (co + [0] * (X + 1))[:X + 1]

def _inv(c, X):
    assert c[0] == 1
    b = [0] * (X + 1); b[0] = 1
    for n in range(1, X + 1):
        b[n] = -sum(c[k] * b[n - k] for k in range(1, min(n, len(c) - 1) + 1))
    return b

def disc_tables(Dmax):
    """A(D) = coefficient of phi_{0,1}/4 and Bm(D) = coefficient of phi_{-2,1} at 4n - j^2 = D,
    from the theta functions, as series in w = q^(1/2)."""
    X = Dmax // 4 + 2; W = 2 * X + 2; R = isqrt(2 * W) + 4
    def sq(kind, j):           # theta_kind(tau, z)^2 at y^j (q^(1/4) removed for kinds 1, 2)
        out = [0] * (W + 1)
        for a in range(-R, R + 1):
            if kind in (1, 2):
                b = j - 1 - a; e = a * a + a + b * b + b
                sgn = -1 if (kind == 1 and (a + b) % 2) else 1
            else:
                b = j - a; e = a * a + b * b
                sgn = -1 if (kind == 4 and (a + b) % 2) else 1
            if e <= W: out[e] += sgn
        return out
    def sq0(kind):             # theta_kind(tau, 0)^2
        one = [0] * (W + 1)
        for a in range(-R, R + 1):
            e = a * a + a if kind == 2 else a * a
            if e <= W: one[e] += (-1 if (kind == 4 and a % 2) else 1)
        return _conv(one, one, W)
    eta6 = [0] * (W + 1); eta6[0] = 1
    for n in range(1, W // 2 + 1):
        for _ in range(6):
            for k in range(W, 2 * n - 1, -1): eta6[k] -= eta6[k - 2 * n]
    ieta6 = _inv(eta6, W)
    A, Bm = {}, {}
    for j in (0, 1):
        ph = _conv(sq(1, j), ieta6, W)
        tot = [Fr(0)] * (W + 1)
        for kind in (2, 3, 4):
            den = sq0(kind); lead = den[0]
            assert all(v % lead == 0 for v in den)
            pr = _conv(sq(kind, j), _inv([v // lead for v in den], W), W)
            for k in range(W + 1): tot[k] += Fr(pr[k], lead)
        for k in range(W + 1):
            if k % 2:
                assert ph[k] == 0 and tot[k] == 0
                continue
            D = 2 * k - j * j
            if -1 <= D <= Dmax: Bm[D] = ph[k]; A[D] = tot[k]
    assert (Bm[-1], Bm[0], Bm[3], Bm[4]) == (1, -2, 8, -12)
    assert (A[-1], A[0], A[3], A[4]) == (Fr(1, 4), Fr(10, 4), -16, 27)
    return A, Bm

class Expo:
    """Exponent c(u, l, r) of (1 - q^l x^u y^r) in Phi~ (x = e(sigma/N), q = e(rho), y = e(v)),
    from the twisted genera F^(0,0) = 8A/N, F^(0,s), F^(r,rk') of Appendix app:computation."""
    def __init__(self, N, Hmax):
        self.N = N
        self.A, self.Bm = disc_tables(4 * Hmax // N + 8)
        c = Fr(24, N - 1)
        def s1(n): return sum(d for d in range(1, n + 1) if n % d == 0)
        e = [Fr(1)] + [c * (s1(n) - (N * s1(n // N) if n % N == 0 else 0)) for n in range(1, Hmax + 3)]
        den = 1
        for v in e: den = den * v.denominator // gcd(den, v.denominator)
        E = [int(v * den) for v in e]
        Dm = 4 * Hmax + 8
        self.BEt = {}
        for rho in (-1, 0, 1, 2):         # sum_t e_t b(D - 4t), D = 4i + rho
            g = [int(self.b(4 * i + rho)) for i in range(Dm // 4 + 2)]
            for i, v in enumerate(_conv(E, g, Dm // 4 + 1)): self.BEt[4 * i + rho] = Fr(v, den)
        self.BENt = {}
        for j in range(isqrt(4 * Hmax // N + 1) + 3):   # sum_s e_{k - N s} b(4s - j^2)
            g = [int(self.b(4 * s - j * j)) for s in range(Hmax // N + 2)]
            for rho in range(N):
                er = E[rho::N]
                for t, v in enumerate(_conv(er, g, len(er) - 1)): self.BENt[(rho + N * t, j)] = Fr(v, den)
        self.cache = {}
    def a(self, D): return self.A.get(D, 0) if D >= -1 else 0
    def b(self, D): return -self.Bm.get(D, 0) if D >= -1 else 0
    def __call__(self, u, l, j):
        key = (u, l, j)
        if key in self.cache: return self.cache[key]
        N = self.N; H = u * l; dl = (l % N == 0)
        if u % N == 0:
            D = 4 * (H // N) - j * j
            F00 = Fr(8, N) * self.a(D)
            F01 = Fr(8, N * (N + 1)) * self.a(D) - Fr(2, N + 1) * (self.BEt[D] if D >= -1 else 0)
            tot = F00 + ((N if dl else 0) - 1) * F01
        else:
            tot = Fr(2, N + 1) * self.BENt[(H, abs(j))]
            if dl: tot += Fr(8, N + 1) * self.a(4 * (H // N) - j * j)
        assert tot.denominator == 1
        self.cache[key] = int(tot)
        return int(tot)

def chamber_coefficients(N, Lq, Ku):
    """c(m, n, l) for -1 <= m < Lq, -1 <= n < Ku: the expansion of
    1/Phi~ = (q x y)^(-1) (1 - y^(-1))^(-2) prod' (1 - q^l x^u y^r)^(-c(u,l,r)) with |y| > 1."""
    ex = Expo(N, Lq * Ku + 2)
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
                    assert abs(r) <= R
                    out[(b, r)] = Fr(int(v.p), int(v.q))
        return out
    Lg = [dict() for _ in range(Lq + 1)]
    sym = 0
    for l in range(Lq + 1):
        for u in range(Ku + 1):
            if l == u == 0: continue
            jm = isqrt(4 * u * l // N + 1) + 1
            for r in range(-jm, jm + 1):
                c = ex(u, l, r)
                assert c == ex(l, u, r), ('q <-> x symmetry', u, l, r)
                sym += 1
                if not c: continue
                t = 1
                while t * l <= Lq and t * u <= Ku:
                    key = (t * u, t * r)
                    Lg[t * l][key] = Lg[t * l].get(key, 0) + Fr(c, t)
                    t += 1
    assert [ex(0, 0, r) for r in range(-4, 0)] == [0, 0, 0, 2]
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
        assert v.denominator == 1
        return int(v)
    return coef, sym

# ---------------------------------------------------------------- B1: orbits
def comp(F, g):
    A, B, C = F; a, b, c, d = g
    return (A*a*a + B*a*b + C*b*b, 2*A*a*c + B*(a*d + b*c) + 2*C*b*d, A*c*c + B*c*d + C*d*d)

def mul(g, h):
    a, b, c, d = g; e, f, x, y = h
    return (a*e + b*x, a*f + b*y, c*e + d*x, c*f + d*y)

def reduce(f):
    """V with f o V = F reduced: -A < B <= A <= C, and B >= 0 if A = C."""
    V = (1, 0, 0, 1); F = f
    while True:
        A, B, C = F
        k = -((B + A - 1) // (2 * A))
        if k: M = (1, 0, k, 1); F = comp(F, M); V = mul(M, V)
        A, B, C = F
        if A > C: M = (0, 1, -1, 0); F = comp(F, M); V = mul(M, V); continue
        if A == C and B < 0: M = (0, 1, -1, 0); F = comp(F, M); V = mul(M, V)
        return F, V

def reduced_forms(delta):
    out = []
    for A in range(1, isqrt(delta // 3) + 1):
        for B in range(-A + 1, A + 1):
            if (delta + B * B) % (4 * A): continue
            C = (delta + B * B) // (4 * A)
            if C > A or (C == A and B >= 0): out.append((A, B, C))
    return out

def automorphs(F):
    return [h for h in ((a, b, c, d) for a in range(-2, 3) for b in range(-2, 3)
            for c in range(-2, 3) for d in range(-2, 3) if a * d - b * c == 1) if comp(F, h) == F]

def row_class(N, row, H):
    c, d = row
    return min(((c * a_ + d * c_) % N, (c * b_ + d * d_) % N) for a_, b_, c_, d_ in H)

def label(N, ch):
    m, n, l = ch
    F, V = reduce((m, N * l, N * n))
    a, b, c, d = V; g = (d, -b, -c, a)
    assert comp(F, g) == (m, N * l, N * n)
    return F, row_class(N, (g[2], g[3]), automorphs(F))

def lift(N, c, d):
    """A matrix in SL_2(Z) with bottom row congruent to (c, d) mod N."""
    for t in range(N * N + 1):
        for cc, dd in ((c, d + N * t), (c + N * t, d)):
            if gcd(cc, dd) == 1:
                if cc == 0: return (dd, 0, 0, dd)
                for a in range(abs(cc) + 1):
                    if (a * dd - 1) % cc == 0: return (a, (a * dd - 1) // cc, cc, dd)
    raise AssertionError((N, c, d))

def orbit_labels(N, D):
    labs = set()
    for F in reduced_forms(N * D):
        H = automorphs(F)
        for c in range(N):
            for d in range(N):
                if (c, d) == (0, 0): continue
                g = lift(N, c, d)
                _, B2, C2 = comp(F, g)
                if B2 % N or C2 % N: continue
                labs.add((F, row_class(N, (c, d), H)))
    return labs

# ---------------------------------------------------------------- B2: the set A
def inA(N, m, n, l):
    D = 4 * m * n - N * l * l
    if not (m > 0 and n > 0 and D > 0): return False
    m1 = N
    while m1 * m1 * D < N * m * m:
        for M1 in (m1, -m1):
            # a violating j has |j m + l M1| < m, so j lies in the open interval of length 2
            # about -l M1 / m and hence in floor(-l M1 / m) + {0, 1}; integer floor, no rounding
            c = (-l * M1) // m
            for j in range(c - 1, c + 3):
                if j % 2 == 0 or ((j * j - 1) // 4) % M1: continue
                if N * (j * m + l * M1) ** 2 + M1 * M1 * D < N * m * m: return False
        m1 += N
    return True

def candidates(N, D, Mb, nb, Lb):
    for m in range(1, Mb + 1):
        for l in range(2 * m):
            if (D + N * l * l) % (4 * m) == 0: yield (m, (D + N * l * l) // (4 * m), l)
    for n in range(1, nb + 1):
        for l in range(Lb):
            if (D + N * l * l) % (4 * n) == 0:
                m = (D + N * l * l) // (4 * n)
                if m > Mb and l < 2 * m: yield (m, n, l)

def check_least(N, Dmax):
    """On a range: in every orbit the charges of least m lie in A."""
    n = 0
    for D in range(1, Dmax + 1):
        labs = orbit_labels(N, D); seen = {}
        m = 0
        while len(seen) < len(labs):
            m += 1
            for l in range(2 * m):
                if (D + N * l * l) % (4 * m) == 0:
                    ch = (m, (D + N * l * l) // (4 * m), l); lb = label(N, ch)
                    if seen.get(lb, (m,))[0] == m:
                        seen[lb] = ch
                        assert inA(N, *ch), ('least charge outside A', N, D, ch)
                        n += 1
    return n

SEARCH = {2: (60, 12, 80), 3: (60, 12, 80), 5: (400, 0, 0), 7: (400, 0, 0)}
FIRST = {2: [360, 2164, 4352, 18944, 36024, 125860, 233080, 705536],
         3: [27, 540, 378, 3294, 2646, 16200, 14094, 68418],
         5: [20, 100, 125, 460, 588, 1720, 2191, 5544],
         7: [3, 36, 18, 138, 81, 444, 273, 1242]}

def run(N, Dmax):
    t0 = time.time()
    Mb, nb, Lb = SEARCH[N]
    reps = {}; norb = 0
    for D in range(1, Dmax + 1):
        labs = orbit_labels(N, D); norb += len(labs)
        found = {}
        for ch in candidates(N, D, Mb, nb, Lb):
            if inA(N, *ch):
                lb = label(N, ch)
                assert lb in labs, ('label outside the enumeration', N, D, ch)
                found.setdefault(lb, []).append(ch)
        missing = [lb for lb in labs if lb not in found]
        assert not missing, ('orbit without a charge in A', N, D, missing[:3])
        for lb, chs in found.items():
            chs.sort(key=lambda c: (min(c[0], c[1]), max(c[0], c[1])))
            reps[(D, lb)] = chs
    Lq = max(min(c[0][0], c[0][1]) for c in reps.values()) + 1
    Ku = max(max(c[0][0], c[0][1]) for c in reps.values()) + 1
    coef, sym = chamber_coefficients(N, Lq, Ku)
    def dstar(ch):
        m, n, l = ch
        a, b = (m, n) if m <= n else (n, m)
        return (-1) ** (l + 1) * coef(a, b, l)
    vals = {}; agree = 0
    for key, chs in reps.items():
        v = dstar(chs[0])
        for ch in chs[1:]:
            if min(ch[0], ch[1]) < Lq and max(ch[0], ch[1]) < Ku:
                assert dstar(ch) == v, ('representatives of one orbit disagree', N, key, chs[0], ch)
                agree += 1
        vals[key] = (chs[0], v)
    assert all(v > 0 for _, v in vals.values()), [k for k, (c, v) in vals.items() if v <= 0][:5]
    # m = 1, Delta = 4n - N l^2 with 0 <= l <= 1: the first values
    first = []
    for D in range(1, Dmax + 1):
        for l in (0, 1):
            if (D + N * l * l) % 4 == 0 and len(first) < 8 and (D + N * l * l) // 4 < Ku:
                n = (D + N * l * l) // 4
                if n > 0: first.append(dstar((1, n, l)))
    k = min(len(first), 8)
    assert first[:k] == FIRST[N][:k], (first, FIRST[N])
    nleast = check_least(N, Dmax)
    least = min(vals.items(), key=lambda kv: kv[1][1])
    h = hashlib.sha256(json.dumps(sorted([[k[0], str(k[1]), list(c), str(v)] for k, (c, v) in vals.items()])).encode()).hexdigest()
    res = {'N': N, 'Delta_max': Dmax, 'orbits': norb, 'box': [Lq - 1, Ku - 1],
           'exponents_checked_symmetric': sym, 'least_charges_in_A_checked': nleast, 'second_representatives_agreeing': agree,
           'first_values_m1': first[:k], 'least': {'Delta': least[0][0], 'charge_m_n_l': list(least[1][0]), 'dstar': least[1][1]},
           'all_positive': True, 'sha256_of_values': h, 'seconds': round(time.time() - t0, 1)}
    print(json.dumps(res), flush=True)
    res['values'] = sorted([[k[0], c[0], c[1], c[2], v] for k, (c, v) in vals.items()])
    return res

if __name__ == '__main__':
    args = sys.argv[1:]; out = None
    if '--out' in args:
        i = args.index('--out'); out = args[i + 1]; args = args[:i] + args[i + 2:]
    vout = None
    if '--values' in args:
        i = args.index('--values'); vout = args[i + 1]; args = args[:i] + args[i + 2:]
    jobs = [tuple(int(x) for x in a.split(':')) for a in args] or [(2, 200), (3, 200), (5, 20), (7, 20)]
    results = [run(N, D) for N, D in jobs]
    if vout: json.dump([{'N': r['N'], 'values': [[d, m, n, l, str(v)] for d, m, n, l, v in r['values']]} for r in results], open(vout, 'w'))
    for r in results: del r['values']
    if out: json.dump(results, open(out, 'w'), indent=1)
    print('all checks passed')
