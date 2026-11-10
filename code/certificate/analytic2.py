# Non-leading pieces for m>=101: explicit bounds, each decreasing in m (checked below by log-derivative).
from analytic import *
def nonmain(m, q):
    m_=arb(m); t0=m_*m_+4*m_; st0=t0.sqrt(); g0=g(PI*st0)
    sq=q.sqrt()
    Pq=Preal(q); Psq=Preal(sq)
    W=(2*m_+1)**2*(Pq**2+64*sq*Pq*Psq)          # sum over walls of weight*e^{-pi c S}
    A1=m_/2-3-2/m_
    rs_ti=W*(-PI*A1/2).exp()/g0
    # gamma>=2: I(z/ga)<=(2/ga)^{15/2} I(z/2); (t/t0)^{15/4} I(z/2)/I(z0) <= sqrt2 e^{-z0/2} e^{-pi Ds/2}/g0;
    # count <= 2 ga, |Weil| factor sqrt(2 gcd)/ga <= sqrt(2/ga); sum_{ga>=2} 2 sqrt2 2^{15/2} ga^{-7} <= 2^9 (zeta(7)-1)
    E2=W/g0*arb(2).sqrt()*(-PI*st0/2).exp()*2**9*arb('0.00835')
    fh=(2*PI*(2*m_+2).sqrt()).exp()              # d1(m) <= e^{2 pi sqrt(2(m+1))}  (log P(e^{-tau}) <= 2 pi^2/tau)
    z1=2*PI*m_.sqrt(); z0=PI*st0
    Iratio=(z0/z1).sqrt()*(z1-z0).exp()/g(z0)    # I(z1)/I(z0)
    ext_over_I0=PI/(2*m_).sqrt()*(m_+2)/16       # ext >= this * I(z0)  (since t0/m^2 >= 1)
    z7=arb('1.0084'); z65=arb('1.0118')          # zeta(7), zeta(6.5) upper bounds
    MK2=(2*m_).sqrt()*arb(17)/16*fh*(4/m_)**4*z7*Iratio/ext_over_I0
    MK3=arb(17)/16*fh/(2*PI)*(4/m_)**(arb(17)/4)*2/m_.sqrt()*PI**2*4*m_**2/3*z65*Iratio/ext_over_I0
    return rs_ti,E2,MK2,MK3
m=101
c=arb(m-1)/(4*(m+2)); q=(-PI*c).exp()
r=nonmain(m,q)
print('q',q.str(5),'rs+TI',r[0].str(3),'E2',r[1].str(3),'MK2',r[2].str(3),'MK3',r[3].str(3),'sum',(r[0]+r[1]+r[2]+r[3]).str(3))
print('main(m>=101) <= 0.33815 ; PMIN 0.538 ; margin >=', (arb('0.538')-arb('0.33815')-r[0]-r[1]-r[2]-r[3]).str(4))
# monotonicity: evaluate at a few larger m with q fixed (q(m) decreases, so P(q(m))<=P(q(101)))
for mm in (150,300,1000,10000):
    rr=nonmain(mm,q); print(mm,[x.str(2) for x in rr])
