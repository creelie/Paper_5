# Interval proof (Arb balls, no monotonicity input) of:
#  (a) phi(z)=z g'(z)/g(z) is decreasing on [3,inf)   (g(z)=sqrt(2 pi z) e^{-z} I_{15/2}(z))
#  (b) phi(z) > 0 on [3,inf), hence g increasing there
#  (c) z*phi(z) <= 29 on [30,inf), hence g'(z) <= 29/z^2 there (as g<=1)
# phi = nu+1/2 - z + z r,  phi' = z-1-2 nu r - z r^2,  r = I_{nu+1}/I_nu.
from flint import arb, ctx
from math import factorial
ctx.prec=120
nu=arb(15)/2
def r(z): return z.bessel_i(nu+1)/z.bessel_i(nu)
def phi(z): return nu+arb(1)/2-z+z*r(z)
def dphi(z): rr=r(z); return z-1-2*nu*rr-z*rr*rr
def d2phi(z):
    rr=r(z); rp=1-(2*nu+1)*rr/z-rr*rr
    return 1-2*nu*rp-rr*rr-2*z*rr*rp
def dphi_box(x,h):   # mean value form: dphi(c) + [-h/2,h/2]*d2phi(box)
    c=arb(x+h/2); Z=arb(x).union(arb(x+h))
    return dphi(c)+arb(0,h/2)*d2phi(Z)
def phi_box(x,h):
    c=arb(x+h/2); Z=arb(x).union(arb(x+h))
    return phi(c)+arb(0,h/2)*dphi(Z)
ok=True; x=3.0; N=0; zmax=60.0
while x<zmax:
    z1=arb(x); est=float(dphi(z1).mid())
    h=min(0.05,0.2*abs(est)); h=max(h,1e-4); h=min(h,zmax-x)
    v=dphi_box(x,h)
    if not v<0:
        h/=4; v=dphi_box(x,h)
        if not v<0: ok=False; print('fail a',x,v); break
    x=min(x+h,zmax); N+=1
print('(a) phi decreasing on [3,60]:',ok,N,'boxes')
# (b) phi(60)>0 (then phi>0 on [3,60] by (a), and phi=h(w)+c>0 for z>=60)
print('(b) phi(60) =',phi(arb(60)).str(6))
# (c) on [30,60]: z*phi(z) <= 28
ok3=True; x=30.0
while x<60:
    v=arb(x+0.05)*phi(arb(x))   # phi decreasing by (a), phi>0
    if not v<29: ok3=False; print('fail c',x,v); break
    x+=0.05
print('(c) z phi <= 29 on [30,60]:',ok3, ' z phi(30)=',(arb(30)*phi(arb(30))).str(6))
# z>=60: g(z)= p(-w) + e^{-2z} p(w), w=1/(2z), p(w)=sum_{k=0}^7 (7+k)!/(k!(7-k)!) w^k  (Watson 3.71 / I_{n+1/2})
a=[arb(factorial(7+k))/(factorial(k)*factorial(7-k)) for k in range(8)]
def p(z,d=0):
    if d==0: return sum((a[k]*z**k for k in range(8)),arb(0))
    if d==1: return sum((k*a[k]*z**(k-1) for k in range(1,8)),arb(0))
    return sum((k*(k-1)*a[k]*z**(k-2) for k in range(2,8)),arb(0))
# main part h(w) = phi_0 = z d/dz log p(-w) = w p'(-w)/p(-w)  (since dw/dz=-w/z)
ok2=True; M=600; minhp=None; maxzphi=None
for j in range(M):
    W=arb(j/(120*M)).union(arb((j+1)/(120*M)))
    hp=(p(-W,1)-W*p(-W,2))/p(-W) + W*p(-W,1)**2/p(-W)**2      # h'(w)
    if not hp>0: ok2=False; print('fail w',j,hp); break
    lo=hp.lower(); minhp=lo if minhp is None else min(minhp,lo)
    zphi=p(-W,1)/(2*p(-W))                                     # z*h(w) = h/(2w)
    up=zphi.upper(); maxzphi=up if maxzphi is None else max(maxzphi,up)
print("h'(w)>0 on (0,1/120]:",ok2,'min',minhp,' max z*h:',maxzphi)
# correction c(z)= z d/dz log(1+eps), eps=e^{-2z}p(w)/p(-w): |c|,|c'| <= 10*z^2*e^{-2z} for z>=60 (crude):
z=arb(60); W=1/(2*z); eps=(-2*z).exp()*p(W)/p(-W)
print('eps(60) <=',eps.str(3),' bound 10 z^2 e^{-2z} at 60:',(10*z*z*(-2*z).exp()).str(3),
      '  vs  min h\'/(2z^2) at z=60:', (arb(minhp)/(2*z*z)).str(4))
