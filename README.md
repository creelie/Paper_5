# Every Single-Centred Black Hole Index of the Z2 CHL Model Is Positive

Paper and verification code by Deep Bhattacharjee.

**Proved.** For every charge (Q, P) of the Z2 CHL model with D = Q²P² − (Q·P)² > 0, the single-centred index d*(Q, P) is positive. The proof combines an exact level-two Rademacher expansion of the finite parts (proved in the paper), interval certificates, and an exact computation for D ≤ 800.

**Layout.** `paper/` — LaTeX source (`main.tex`, `parts/`), figures (`paper/figures/`) and PDF. `code/chl_coefficients.py` — exact orbit computation for D ≤ 800; `code/verify_table.py` — the table of values. `code/certificate/` — the ball-arithmetic certificates. `code/proof_checks/` — independent checks of the proof of the expansion. Outputs (`*.out`) are stored next to the programs.

**Running the checks** (Python 3 with `python-flint`, `mpmath`, `numpy`):

```
python3 code/chl_coefficients.py 2:800    # exact indices of all orbits with D <= 800
python3 code/verify_table.py              # Table of values
python3 code/certificate/cert.py -x 4 1   # certificate at Q^2 = 1, x0 = 4 (cert.py 8 9 ... for x0 = 1)
python3 code/proof_checks/polar_check.py  # likewise polar0_check.py, transform_check.py, chk_*.py, hp.py
python3 code/proof_checks/expansion_check.py 1 8 0 4 1 1   # truncated Rademacher series at Q^2 = 1
cd paper && latexmk -pdf main.tex
```

DOI (all versions): [10.5281/zenodo.24094931](https://doi.org/10.5281/zenodo.24094931)
