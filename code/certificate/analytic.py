# Bounds for m >= 101 (all in Arb).
import sys
from flint import arb, ctx
ctx.prec=96
import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import polar
D1=polar.D1
PI=arb.pi()
NU=arb(15)/2
SMAX=300
p2=[0]*(SMAX+2)
for S in range(0,SMAX+2):
    p2[S]=sum(D1(a)*D1(S-2-a) for a in range(-1,S))
def g(y): return (2*PI*y).sqrt()*(-y).exp()*y.bessel_i(NU)
def Preal(q,N=400):
    v=arb(1)
    for n in range(1,N+1): v*=((1-q**n)*(1-q**(2*n)))**8
    return (16*q**(N+1)/((1-q)*(1-q**(N+1)))).exp()/v
def phi(S,m):
    m=arb(m); S_=arb(S)
    etam=S_*(S_+4)/(2*m+4)                 # S - Delta_s <= etam
    # Delta_s - S <= etap: for S <= (m+4)/2, t >= m(m+4-2S) since 4ab <= (a+b)^2 = (S-2)^2;
    # otherwise only Delta_s <= sqrt(t0) <= m+2 is used.
    if 2*S<=m+4:
        etap=((m*(m+4)).sqrt()-(m*(m+4-2*S_)).sqrt()-S_).max(arb(0))
    else:
        etap=m+2-S_
    ytl2=m*m+4*m-2*m*S_
    eps=arb(1)
    if ytl2>0:
        ytl=PI*ytl2.sqrt()
        if ytl>30:
            epsg=PI*(S_+etap)*29/(ytl**2*g(ytl))
            eps=(8*S_/m+epsg).min(arb(1))
    return eps+(PI*etam).exp()-1+(PI*etap).min(arb(1))
def main_bound(m,prt=False):
    u=(-PI).exp(); tot=arb(0)
    for S in range(1,min(m+1,SMAX)+1):
        tot+=arb(p2[S])*u**S*phi(S,m)
    # missing pairs (S >= m/2+2) and S > m+1: model terms
    for S in range(max(1,(m+3)//2+1),SMAX+1):
        tot+=arb(p2[S])*u**S
    q=(-arb('0.4')*PI).exp(); rho=(-arb('0.485')*PI).exp()/q
    tot+=2*Preal(q)**2*rho**(SMAX+1)/(1-rho)
    return tot
if __name__=="__main__":
    worst=None
    for m in range(101,301):
        b=main_bound(m)
        if worst is None or b>worst[0]: worst=(b,m)
        if m%20==1: print(m,b.str(5,radius=False),flush=True)
    print('max over 101..300:',worst[0].str(5),'at',worst[1])
    # m>=300: all S<=300 terms decreasing in m -> bound at m=300 with S range full (no missing-sum change since S<=300<=m)
    print('m>=300 bound:',main_bound(300).str(5))
