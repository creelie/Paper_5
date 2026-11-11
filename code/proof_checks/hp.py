import mpmath as mp, sys
mp.mp.dps=35
def theta_direct(t,m,ga,g):
    # sum_p exp(-pi (g+2m ga p)^2 t/(2m)), enough terms
    a=2*m*ga; P=int(mp.sqrt(60/(mp.pi*a*a*t/(2*m))))+3
    return mp.fsum(mp.exp(-mp.pi*(g+a*p)**2*t/(2*m)) for p in range(-P-abs(g)//a-1,P+abs(g)//a+2))
def theta_dual(t,m,ga,g):
    # Poisson: = (1/sqrt(A)) sum_k e(k g/(2m ga)) exp(-pi k^2/A), A=2m ga^2 t
    A=2*m*ga**2*t; c=mp.mpf(g)/(2*m*ga)
    K=int(mp.sqrt(60*A/mp.pi))+3
    return mp.re(mp.fsum(mp.expjpi(2*k*c)*mp.exp(-mp.pi*k*k/A) for k in range(-K,K+1)))/mp.sqrt(A)
def theta(t,m,ga,g):
    return theta_dual(t,m,ga,g) if 2*m*ga**2*t<1 else theta_direct(t,m,ga,g)
def f(m,ga,g,u):
    if g%(2*m*ga)==0:
        if abs(u)<mp.mpf(10)**-12: return -mp.pi**2/3
        return mp.pi**2/mp.sinh(mp.pi*u/ga)**2-mp.mpf(ga)**2/u**2
    return mp.pi**2/mp.sinh(mp.pi*u/ga-1j*mp.pi*mp.mpf(g)/(2*m*ga))**2
def lhs(m,ga,g,X):
    # t=s^2 removes t^{-1/2}
    return mp.quad(lambda s: 2*s*theta(s*s,m,ga,g)*(s*s+X)**mp.mpf(-1.5),[0,0.05,0.2,0.5,1,2,4,mp.inf])
def rhs(m,ga,g,X):
    a=2/mp.sqrt(X) if g%(2*m*ga)==0 else 0
    return a-1/(mp.sqrt(2*m)*mp.pi*ga**2*X)*mp.quad(lambda u: mp.exp(-2*mp.pi*m*X*u*u)*f(m,ga,g,u),[-mp.inf,-2,-1,-0.5,0,0.5,1,2,mp.inf])
if __name__=="__main__":
 for (m,ga,g,X) in [(1,1,1,mp.mpc(0.7,0.3)),(1,2,0,mp.mpc(0.5,-0.4)),(2,3,5,mp.mpc(1.3,0.9)),(3,2,-6,mp.mpc(0.4,0.2))]:
     L=lhs(m,ga,g,X); R=rhs(m,ga,g,X); print(m,ga,g,X,mp.nstr(L,20),mp.nstr(abs(L-R),3),flush=True)
