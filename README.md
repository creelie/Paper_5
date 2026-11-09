# Single-Centred Dyon Indices of the Z2 CHL Model

Paper and verification code for Sen's positivity conjecture in the Z2 CHL model of N = 4 string theory, by Deep Bhattacharjee.

- `paper/` — LaTeX source, figures and PDF; `TexmdpiCHLPositivity.zip` is the TeX package.
- `code/` — exact-arithmetic checks (Python 3 with `python-flint`):
  `python3 code/verify_table.py` reproduces Table 1 and checks positivity on a box of charges;
  `python3 code/chl_coefficients.py 2:200` checks every duality orbit with discriminant up to 200.

Proved: positivity for every charge with Q² = 1 or 2 and its duality orbit. Open: Q² ≥ 3.

Contact for GitHub matters: egspl@proton.me.
