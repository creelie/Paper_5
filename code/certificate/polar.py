# polar data of the Z2 CHL finite parts at the two cusps, from the wall formulas (P0), (Pinf)
from fractions import Fraction as Fr
from math import gcd
import functools
JMAX=3000
_d1=None
def d1tab(J=JMAX):
    global _d1
    if _d1 is None or len(_d1)<J+1:
        d=[0]*(J+1); d[0]=1
        for r in range(1,J+1):
            e=16 if r%2==0 else 8
            for _ in range(e):
                for k in range(r,J+1): d[k]+=d[k-r]
        _d1=d
    return _d1
def D1(x):  # coefficient of q^x in 1/(eta^8 eta^8(2.)), x>=-1 integer
    if x<-1: return 0
    return d1tab()[x+1]
def D2(t):  # d1(2t), t in Z/2
    u=2*t
    if u.denominator!=1: return 0
    return D1(int(u))
def walls(m,nu,lt,R=None):
    """matrices G=(p q;r s), r,s>=1, ps-qr=1, q/s <= lt/2m < p/r ; returns (r,s,mG,nG,lG)"""
    x0=Fr(lt,2*m)
    if R is None: R=2*m+6
    out=[]
    for s in range(1,R+1):
        for r in range(1,R+1):
            if gcd(r,s)!=1: continue
            # q = -r^{-1} mod s ; q/s <= x0 < q/s+1/(rs)
            if s==1: q0=0
            else: q0=(-pow(r,-1,s))%s
            # choose q = q0 + s*t with q/s <= x0 < q/s + 1/(rs)
            import math
            t=math.floor((x0*s-q0)/s)
            q=q0+s*t
            if not (Fr(q,s)<=x0<Fr(q,s)+Fr(1,r*s)): continue
            p=(1+q*r)//s
            assert p*s-q*r==1
            nG=s*s*nu+q*q*m-q*s*lt; mG=r*r*nu+p*p*m-p*r*lt; lG=-2*r*s*nu-2*p*q*m+(p*s+q*r)*lt
            out.append((r,s,mG,nG,lG))
    return out
def polarS(m):
    """dict (Dt, lt) -> 16*coefficient, 0<=lt<=m, Dt=4m*n-lt^2<0"""
    res={}
    for lt in range(0,m+1):
        n=-1
        while 4*m*n-lt*lt<0:
            v=0
            for (r,s,mG,nG,lG) in walls(m,Fr(n),lt):
                if mG>=-1 and nG>=-1: v+=lG*D1(int(mG))*D1(int(nG))
            if v: res[(4*m*n-lt*lt,lt)]=int(v)
            n+=1
    return res
def polarI(m):
    """dict (Dt, lt) -> coefficient of phi^F_m at infinity, Dt=2m e-lt^2<0 (nu=e/2)"""
    res={}
    for lt in range(0,m+1):
        e=-1
        while 2*m*e-lt*lt<0:
            nu=Fr(e,2); v=0
            for (r,s,mG,nG,lG) in walls(m,nu,lt):
                if r%2==0 and mG.denominator==1 and mG>=-1 and 2*nG>=-1: v+=lG*D1(int(mG))*D2(nG)
                if s%2==0 and nG.denominator==1 and nG>=-1 and 2*mG>=-1: v+=lG*D2(mG)*D1(int(nG))
            if v: res[(2*m*e-lt*lt,lt,e%2)]=int(v)
            e+=1
    return res
if __name__=="__main__":
    import sys
    m=int(sys.argv[1])
    print(sorted(polarS(m).items())); print(sorted(polarI(m).items()))
