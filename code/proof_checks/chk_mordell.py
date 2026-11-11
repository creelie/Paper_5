import mpmath as mp
mp.mp.dps=25
def theta(t,m,ga,g,N=60):
    return mp.fsum(mp.e**(-mp.pi*(g+2*m*ga*n)**2*t/(2*m)) for n in range(-N,N+1))
def f(m,ga,g,u):
    if g%(2*m*ga)==0:
        if abs(u)<1e-8: return -mp.pi**2/3
        return mp.pi**2/mp.sinh(mp.pi*u/ga)**2-mp.mpf(ga)**2/u**2
    return mp.pi**2/mp.sinh(mp.pi*u/ga-1j*mp.pi*g/(2*m*ga))**2
for (m,ga,g,X) in [(1,1,1,mp.mpc(0.7,0.4)),(2,3,5,mp.mpc(0.3,-0.8)),(1,2,0,mp.mpc(1.1,0.5)),(3,1,-2,mp.mpc(0.5,0.2))]:
    lhs=mp.quad(lambda t: theta(t,m,ga,g)*(t+X)**(-1.5),[0,0.01,1,mp.inf])
    rhs=(2*X**-0.5 if g%(2*m*ga)==0 else 0) - 1/(mp.sqrt(2*m)*mp.pi*ga**2*X)*mp.quad(lambda u: mp.e**(-2*mp.pi*m*X*u**2)*f(m,ga,g,u),[-mp.inf,0,mp.inf])
    print(m,ga,g,mp.nstr(lhs,15),mp.nstr(rhs,15))
# Bessel line integral
A,B,k,r=mp.mpf(0.3),mp.mpf(0.7),2,mp.mpf(6.5); c=mp.mpf(0.4)
li=mp.quad(lambda s: (c+1j*s)**(-r-2)*mp.e**((2*mp.pi/k)*(A*(c+1j*s)+B/(c+1j*s))),[-mp.inf,-10,0,10,mp.inf])/(2*mp.pi)
print(li, (A/B)**((r+1)/2)*mp.besseli(r+1,4*mp.pi*mp.sqrt(A*B)/k))
