"""Write parts/margins_table.tex from the certificate outputs."""
import os, re
here = os.path.dirname(os.path.abspath(__file__))
cdir = os.path.join(here, '..', '..', 'code', 'certificate')
rows = {}
for fn in ('cert_m3_7.out', 'cert_m8_100.out'):
    for line in open(os.path.join(cdir, fn)):
        m = re.match(r'(\d+) .*total=(\S+) margin (\S+)', line)
        if m:
            rows[int(m.group(1))] = (float(m.group(2)), float(m.group(3)))
ms = sorted(rows)
def fmt(e):
    if e < 1e-3:
        a, b = ('%.2e' % e).split('e')
        return r'$%s\cdot10^{%d}$' % (a, int(b))
    return '%.4f' % e
G = 4
per = -(-len(ms) // G)
out = [r'\begin{table}[ht]', r'\centering', r'\footnotesize', r'\setlength{\tabcolsep}{3.5pt}',
       r'\caption{The certified total error $\epsilon(m,x_0(m))$ of Theorem~\ref{thm:cert} and the margin $0.53816-\epsilon$, for $m=Q^2$.}\label{tab:margins}',
       r'\begin{tabular}{' + '|'.join(['rll'] * G) + '}', r'\hline',
       ' & '.join([r'$m$ & $\epsilon$ & margin'] * G) + r' \\', r'\hline']
for i in range(per):
    cells = []
    for g in range(G):
        k = g * per + i
        if k < len(ms):
            m = ms[k]; e, mg = rows[m]
            cells.append('%d & %s & %.4f' % (m, fmt(e), mg))
        else:
            cells.append(' & & ')
    out.append(' & '.join(cells) + r' \\')
out += [r'\hline', r'\end{tabular}', r'\end{table}']
open(os.path.join(here, '..', 'parts', 'margins_table.tex'), 'w').write('\n'.join(out) + '\n')
print(len(ms), 'rows; max eps', max(v[0] for v in rows.values()), 'min margin', min(v[1] for v in rows.values()))
