#!/usr/bin/env python3
"""Checks Table 1 and Theorem A of the paper on a box of charges of the Z2 CHL model.

Charges: m = Q^2, n = P^2/2, l = Q.P, D = Q^2 P^2 - (Q.P)^2 = 2mn - l^2.
The coefficient c(m, n, l) of 1/Phi~_6 in the chamber R is computed from the product formula in
exact arithmetic (chl_coefficients.chamber_coefficients); on the set A, d* = (-1)^(l+1) c(m, n, l).
Usage: python3 verify_table.py [Kmax]   (default 40)
"""
import sys
from chl_coefficients import chamber_coefficients

K = int(sys.argv[1]) if len(sys.argv) > 1 else 40
coef, _ = chamber_coefficients(2, 4, K + 1)

def dstar(m, n, l):
    return (-1) ** (l + 1) * coef(m, n, l)

TABLE_Q1 = [360, 2164, 4352, 18944, 36024, 125860, 233080, 705536]
TABLE_Q2 = {1: [4352, 67008, 671744, 5189696, 33581312, 190095232, 968248576, 4522209472],
            0: [18944, 198144, 1580672, 10586880, 62035968, 327011328, 1578794128, 7075522048],
            2: [6912, 101376, 992256, 7550400, 48253952, 270195712, 1362534400, 6304638848]}

# Q^2 = 1: D = 2n - l^2, l = D mod 2
got = [dstar(1, (D + D % 2) // 2, D % 2) for D in range(1, 9)]
assert got == TABLE_Q1, got
# Q^2 = 2: D = 4n - l^2 with l = 1, 0, 2
for l, row in TABLE_Q2.items():
    Ds = [4 * j - l * l + (4 if l == 2 else 0) for j in range(1, 9)] if l != 1 else [4 * j - 1 for j in range(1, 9)]
    got = [dstar(2, (D + l * l) // 4, l) for D in Ds]
    assert got == row, (l, got)

count = 0
for m in (1, 2):
    for n in range(1, K + 1):
        for l in range(0, 2 * m):
            if 2 * m * n - l * l > 0:
                assert dstar(m, n, l) > 0, (m, n, l)
                count += 1
print(f'Table 1 reproduced; d* > 0 on all {count} charges with Q^2 in (1, 2), 1 <= P^2/2 <= {K}, 0 <= Q.P < 2Q^2, D > 0')
