# rigorous lower bound for min_{|z|<=e^{-pi}} |P(z)|^2, P(z)=prod (1-z^n)^-8 (1-z^{2n})^-8
from flint import arb, acb, ctx
ctx.prec=80
PI=arb.pi()
N=60
def Pc(z):  # complex, with tail enclosure
    v=acb(1)
    for n in range(1,N+1):
        v*=((1-z**n)*(1-z**(2*n)))**8
    v=1/v
    r=abs(z)
    tail=16*r**(N+1)/((1-r)*(1-r**(N+1)))
    # |log of omitted factors| <= tail  =>  |P - P_N| <= |P_N|(e^tail-1)
    err=abs(v)*(tail.exp()-1)
    return v, err
r0=(-PI).exp()
Pr,_=Pc(acb(r0)); Pr=Pr.real
# log-derivative at r0 (positive terms), P'(r0)
ld=sum((8*n*r0**(n-1)/(1-r0**n)+16*n*r0**(2*n-1)/(1-r0**(2*n)) for n in range(1,400)),arb(0))
ld+=arb(1e-100)
Lip=2*Pr*r0*Pr*ld*arb('1.0001')   # |d|P|^2/dtheta| <= 2|P||zP'| <= 2 P(r0) r0 P'(r0)
M=4000
best=None
for k in range(M+1):
    th=PI*k/M
    z=acb(r0*th.cos(),r0*th.sin())
    v,e=Pc(z)
    lo=(abs(v)-e)**2
    lo=arb(lo.lower())
    if best is None or lo<best: best=lo
bound=best-Lip*PI/(2*M)
print('P(r0)=',Pr,'Lip',Lip,'grid min lower',best,'certified min |P|^2 >=',bound.lower())
