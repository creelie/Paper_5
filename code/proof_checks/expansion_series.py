# The principal-part series M of the Rademacher expansion (eq:radmain) for the Z2 CHL finite part phi^F_m,
# summed over the frames I, T, S of Gamma^0(2); imported by expansion_check.py.
import sys, os, mpmath as mp
from math import gcd
from fractions import Fraction as Fr
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import chl_coefficients as C
mp.mp.dps=30
m=int(sys.argv[1]); GMAX=int(sys.argv[2])
coef,_=C.chamber_coefficients(2,m+2,60)
# d1
J=60; d1=[0]*(J+1); d1[0]=1
for r in range(1,J+1):
    e=16 if r%2==0 else 8
    for _ in range(e):
        for k in range(r,J+1): d1[k]+=d1[k-r]
D1=lambda x: d1[x+1] if -1<=x<J-1 else 0
def wall(m,n,l):
    tot=0
    for r in range(1,4*m+5):
        for s in range(1,4*m+5):
            for q in range(-2*s,2*s+1):
                if (1+q*r)%s: continue
                p=(1+q*r)//s
                x=Fr(l,2*m)-Fr(q,s)
                if not (0<=x<Fr(1,r*s)): continue
                nG=s*s*n+q*q*m-q*s*l; mG=r*r*n+p*p*m-p*r*l; lG=-2*r*s*n-2*p*q*m+(p*s+q*r)*l
                tot+=lG*D1(mG)*D1(nG)
    return tot
# polar data: frame -> list of (nu (Fraction), ltilde mod 2m, coefficient)
pol={'I':[], 'S':[]}
for l in range(2*m):
    for e in range(-1,m+2):
        if 2*m*e-l*l<0:
            lc=min(l,2*m-l); ec=e if lc==l else e-2*l+2*m
            v=coef(m,ec,lc) if ec>=-1 else 0
            if v: pol['I'].append((Fr(e,2),l,mp.mpf(v)))
    for n in range(-1,m+2):
        if 4*m*n-l*l<0:
            lc=min(l,2*m-l); nc=n if lc==l else n-l+m
            v=wall(m,nc,lc)
            if v: pol['S'].append((Fr(n),l,mp.mpf(v)/16))
print('polar I',[(str(a),b,int(c)) for a,b,c in pol['I']]); print('polar S',[(str(a),b,float(c)) for a,b,c in pol['S']])
def frame(G):  # coset of G in Gamma^0(2)\SL2 : reps I, T, S ; G = g R, g has b even
    a,b,c,d=G
    for name,R in (('I',(1,0,0,1)),('T',(1,1,0,1)),('S',(0,-1,1,0))):
        # g = G R^{-1}
        ra,rb,rc,rd=R; Ri=(rd,-rb,-rc,ra)
        gb=a*Ri[1]+b*Ri[3]
        if gb%2==0: return name
def weil(al,ga,de,lt,l):   # psi(Gamma)_{lt,l} from (2.41)
    s=0
    for T in range(ga):
        x=lt-2*m*T
        s+=mp.expjpi(2*(mp.mpf(al)/ga*x*x/(4*m) - mp.mpf(l*x)/(2*m*ga) + mp.mpf(de)/ga*l*l/(4*m)))
    return s/mp.sqrt(2*m*ga*1j)
def egcd_inv(a,n):
    return pow(a,-1,n)
def rad(n,l):
    D=mp.mpf(2*m*n-l*l)   # = 4m nu - l^2, nu=n/2
    nu=mp.mpf(n)/2
    tot=0
    for ga in range(1,GMAX+1):
        for md in range(0,2*ga):   # -delta in [0,2ga)
            de=-md
            if gcd(de,ga)!=1: continue
            # alpha with alpha*delta = 1 mod gamma; beta
            al=pow(de% ga,-1,ga) if ga>1 else 0
            be=(al*de-1)//ga
            G=(al,be,ga,de); Gi=(de,-be,-ga,al)
            fr=frame(Gi)
            data = pol['S'] if fr=='S' else pol['I']
            for (nt,lt,a) in data:
                if fr=='T': a=a*(-1)**int(2*nt)   # psi|T coefficient phase e(nu)
                Dt=4*m*mp.mpf(nt.numerator)/nt.denominator - lt*lt
                ph=mp.expjpi(2*(al*Dt/(4*m)+de*(D/(4*m)))/ga)
                tot+=a*ph*weil(al,ga,de,lt,l)*(2*mp.pi/ga)*(-Dt/D)**(mp.mpf(15)/4)*mp.besseli(mp.mpf(15)/2,mp.pi*mp.sqrt(-Dt*D)/(m*ga))
    tot=tot/2*(-1)**l*mp.sqrt(1j)
    return tot
for n in (4,8,16,30):
    for l in range(0,m+1):
        d=(-1)**(l+1)*coef(m,n,l)
        r=rad(n,l)
        print(n,l,'d*',d,'rad',mp.nstr(r,12),'diff',mp.nstr(d-r,8))
