"""Write a single-file main.tex (all \\input{parts/...} inlined) to the given path."""
import os, re, sys
here = os.path.dirname(os.path.abspath(__file__))
root = os.path.join(here, '..')
def expand(path):
    s = open(path).read()
    return re.sub(r'\\input\{(parts/[^}]*)\}', lambda m: expand(os.path.join(root, m.group(1) + ('' if m.group(1).endswith('.tex') else '.tex'))).rstrip('\n'), s)
open(sys.argv[1], 'w').write(expand(os.path.join(root, 'main.tex')))
