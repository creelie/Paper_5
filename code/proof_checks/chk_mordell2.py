import mpmath as mp
mp.mp.dps=20
def theta(t,m,ga,g):
    # Poisson for small t
    Aa=2*m*ga**2*t; c=mp.mpf(g)/(2*m*ga)
    if t<0.05:
        return mp.re(mp.fsum(mp.expjpi(2*k*c)*mp.e**(-mp.pi*k*k/Aa) for k in range(-6,7))/mp.sqrt(Aa))
    return mp.fsum(mp.e**(-mp.pi*(g+2*m*ga*n)**2*t/(2*m)) for n in range(-40,41))
def f(m,ga,g,u):
    if g%(2*m*ga)==0:
        if abs(u)<1e-6: return -mp.pi**2/3
        return mp.pi**2/mp.sinh(mp.pi*u/ga)**2-mp.mpf(ga)**2/u**2
    return mp.pi**2/mp.sinh(mp.pi*u/ga-1j*mp.pi*g/(2*m*ga))**2
for (m,ga,g,X) in [(1,1,1,mp.mpc(0.7,0.4)),(2,3,5,mp.mpc(0.3,-0.8)),(1,2,0,mp.mpc(1.1,0.5)),(3,1,-2,mp.mpc(0.5,0.2))]:
    lhs=mp.quad(lambda s: 2*s*theta(s*s,m,ga,g)*(s*s+X)**(-1.5),[0,0.1,0.3,1,3,mp.inf])
    rhs=(2*X**-0.5 if g%(2*m*ga)==0 else 0) - 1/(mp.sqrt(2*m)*mp.pi*ga**2*X)*mp.quad(lambda u: mp.e**(-2*mp.pi*m*X*u**2)*f(m,ga,g,u),[-mp.inf,-1,0,1,mp.inf])
    print(m,ga,g,mp.nstr(lhs,15),mp.nstr(rhs,15))
