import mpmath as mp
from math import gcd
mp.mp.dps=25
def th(m,l,u,R=60):
    return mp.fsum(mp.expjpi(2*(l+2*m*n)**2*u/(4*m)) for n in range(-R,R+1))
def W(m,ga,de,al,lt,l):
    s=0
    for T in range(ga):
        x=lt-2*m*T
        s+=mp.expjpi(2*mp.mpf(al*x*x-2*l*x+de*l*l)/(4*m*ga))
    return s/mp.sqrt(2*m*ga*1j)
worst=0
for m in (1,2,3):
  for ga in (1,2,3,4,5):
    for md in range(2*ga):
      de=-md
      if gcd(de,ga)!=1: continue
      al=pow(de%ga,-1,ga) if ga>1 else 0
      be=(al*de-1)//ga
      u=mp.mpc(-al/ga+0.13, 0.9/ga)   # generic point
      Nu=(de*u+be)/(ga*u+al)
      for l in range(2*m):
        lhs=th(m,l,Nu)
        rhs=mp.sqrt(ga*u+al)*mp.fsum(W(m,ga,de,al,lp,l)*th(m,lp,u) for lp in range(2*m))
        worst=max(worst,abs(lhs-rhs)/abs(lhs))
print('worst rel err',worst)
