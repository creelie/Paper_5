# Rigorous (Arb) lower bound for d*/ext at Q^2=m, uniformly in x>=1 and l, under (H),(P0),(Pinf).
import sys, math
from math import gcd
from fractions import Fraction as Fr
from flint import arb, ctx
ctx.prec=96
import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import polar
D1=polar.D1
PI=arb.pi()
NU=arb(15)/2
def I(nu,z): return z.bessel_i(nu)
def g(y): return (2*PI*y).sqrt()*(-y).exp()*I(NU,y)
X0=arb(1)
def Rt(t,t0):   # B(t)/B(t0) at x=X0, t arb
    return (t/t0)**(arb(15)/4)*I(NU,PI*X0*t.sqrt())/I(NU,PI*X0*t0.sqrt())
def amax(*xs):
    v=xs[0]
    for y in xs[1:]: v=v.max(y)
    return v
# ---------- main family (r=s=1) ----------
def main_terms(m):
    out=[]
    for a in range(-1,m+2):
        for b in range(-1,m+2-a):
            l=m-a-b
            if l<1: continue
            t=l*l-4*a*b
            if t<=0 or not (-m<=b-a<m): continue
            out.append((a,b,l,t))
    return out
def hsup(C,al,be):   # sup_{x>=X0} |C e^{-al x}-e^{-be x}|, C,al,be>0 arb; rigorous upper bound
    f1=abs(C*(-al*X0).exp()-(-be*X0).exp())
    # interior critical point of C e^{-al x}-e^{-be x}: value there bounded by max(C e^{-al x*}, e^{-be x*}) <= max(C e^{-al},e^{-be}) if x*>=1
    # so use exact x* when well separated, else crude
    try:
        if (al-be).abs_lower()>1e-12:
            xs=((be/(C*al)).log())/(be-al)
            if xs>X0:
                return f1.max(abs(C*(-al*xs).exp()-(-be*xs).exp()))
            if xs<X0: return f1
    except Exception: pass
    return f1.max(C*(-al*X0).exp()).max((-be*X0).exp())
def err_main(m,t0):
    tot=arb(0)
    st0=t0.sqrt(); g0=g(PI*X0*st0)
    for a,b,l,t in main_terms(m):
        S=a+b+2
        if S==0: continue   # the extreme term itself: finite term = model term = 1
        ta=arb(t); ds=st0-ta.sqrt()
        c0=arb(l)/(m+2)*(ta/t0)**(arb(7)/2)
        gh=g(PI*X0*ta.sqrt())/g0
        al=PI*ds; be=PI*S
        v=hsup(c0,al,be).max(hsup(c0*gh,al,be))
        tot+=arb(D1(a))*arb(D1(b))*v
    return tot
def miss(m):
    S=set((a,b) for a,b,l,t in main_terms(m))
    u=(-PI*X0).exp(); tot=arb(0)
    for a in range(-1,4*m+10):
        for b in range(-1,4*m+10-a):
            if (a,b) not in S: tot+=arb(D1(a))*arb(D1(b))*u**(a+b+2)
    # tail S>K=4m+10: p2(S) <= P(q1)^2 q1^{-S}, q1=sqrt(u)
    K=4*m+10; q1=u.sqrt()
    tot+=Preal(q1)**2*(u/q1)**(K+1)/(1-u/q1)
    return tot
def Preal(q,N=200):
    v=arb(1)
    for n in range(1,N+1): v*=((1-q**n)*(1-q**(2*n)))**8
    tail=16*q**(N+1)/((1-q)*(1-q**(N+1)))
    return tail.exp()/v
# ---------- walls with rs>=2 (S frame) and parity walls (I/T frames) ----------
def delta_hi(m,r,s,a,b):
    """upper bound for delta^2 solving sqrt(d2+A)+sqrt(d2+B)=1/(rs), or None if no solution"""
    A=a/(m*r*r); B=b/(m*s*s); tg=1.0/(r*s)
    lo=max(0.0,-A,-B)
    f=lambda d: math.sqrt(max(d+A,0))+math.sqrt(max(d+B,0))-tg
    lo_a=arb(lo)
    fa=lambda d: (d+A_).max(arb(0)).sqrt()+(d+B_).max(arb(0)).sqrt()-tg_
    A_=arb(a)/(m*r*r); B_=arb(b)/(m*s*s); tg_=arb(1)/(r*s)
    # no solution if f>0 at the least admissible delta^2 (f is increasing); certified in Arb
    if fa(lo_a.max(-A_).max(-B_))>0: return None
    hi=lo+tg*tg+1.0
    for _ in range(200):
        mid=(lo+hi)/2
        if f(mid)>0: hi=mid
        else: lo=mid
    hi=hi*(1+1e-9)+1e-15
    # certify in Arb that hi lies above the root
    assert fa(arb(hi))>0
    return hi
def wall_list(m,parity=False):
    """list of (coef_bound, t_hi): coef_bound >= l_G * weight(a,b)"""
    out=[]; R=2*m+1
    for r in range(1,R+1):
        for s in range(1,R+1):
            if gcd(r,s)!=1: continue
            if not parity and r*s<2: continue
            if parity and r%2 and s%2: continue
            amax_=m/(s*s); bmax_=m/(r*r)
            if parity:
                As=[Fr(k,2) for k in range(-2,int(2*amax_)+1)]; Bs=[Fr(k,2) for k in range(-2,int(2*bmax_)+1)]
            else:
                As=list(range(-1,int(amax_)+1)); Bs=list(range(-1,int(bmax_)+1))
            for a in As:
                for b in Bs:
                    if parity:
                        w=0
                        if r%2==0 and Fr(a).denominator==1 and a>=-1 and 2*b>=-1: w=max(w,D1(int(a))*D1(int(2*b)))
                        if s%2==0 and Fr(b).denominator==1 and b>=-1 and 2*a>=-1: w=max(w,D1(int(2*a))*D1(int(b)))
                        if w==0: continue
                    else:
                        w=D1(a)*D1(b)
                    d2=delta_hi(m,r,s,float(a),float(b))
                    if d2 is None: continue
                    lhi=m-float(a)-float(b)
                    if lhi<=0: continue
                    t_hi=4*m*m*d2
                    out.append((arb(lhi)*arb(w),arb(t_hi)))
    return out
def sum_R(lst,t0,m): return sum((c*Rt(t,t0) for c,t in lst),arb(0))/(m+2)
# ---------- frames ----------
def frame_counts(ga):
    cnt={'I':0,'T':0,'S':0}
    for md in range(0,2*ga):
        de=-md
        if gcd(de,ga)!=1: continue
        al=pow(de%ga,-1,ga) if ga>1 else 0
        be=(al*de-1)//ga
        a,b,c,d=(de,-be,-ga,al)
        for name,R in (('I',(1,0,0,1)),('T',(1,1,0,1)),('S',(0,-1,1,0))):
            ra,rb,rc,rd=R; Ri=(rd,-rb,-rc,ra)
            if (a*Ri[1]+b*Ri[3])%2==0: cnt[name]+=1; break
    return cnt
def certify(m,G=12,verbose=True):
    t0=arb(m*m+4*m); st0=t0.sqrt()
    MT=main_terms(m)
    Smain=[(arb(l)*arb(D1(a))*arb(D1(b)),arb(t)) for a,b,l,t in MT]
    RS=wall_list(m,False); IT=wall_list(m,True)
    E_main=err_main(m,t0); E_miss=miss(m)
    E_rs=sum_R(RS,t0,m)
    E_TI=2*16*sum_R(IT,t0,m)
    # gamma>=2
    I0=I(NU,PI*X0*st0)
    def sumBg(lst,ga): return sum((c*(t/t0)**(arb(15)/4)*I(NU,PI*X0*t.sqrt()/ga) for c,t in lst),arb(0))/(I0*(m+2))
    SS=Smain+RS
    E2=arb(0)
    for ga in range(2,G+1):
        c=frame_counts(ga); gg=gcd(m,ga)
        E2+=(2*arb(gg)).sqrt()/ga*(c['S']*sumBg(SS,ga)+(c['I']+c['T'])*16*sumBg(IT,ga))
    X1=sum_R(SS,t0,m); Y1=16*sum_R(IT,t0,m)
    E2+=2*arb(2).sqrt()*(X1+Y1)/(6*arb(G)**6)
    # mock lines (relative to ext at x=1, D=m^2)
    D=arb(m*m)*X0*X0; fh=arb(D1(m))
    ext=PI/(2*arb(m)).sqrt()*arb(m+2)/16*(t0/D)**(arb(15)/4)*I0
    MK2=arb(0); MK3=arb(0)
    for ga in range(1,G+1):
        c=frame_counts(ga); gg=gcd(m,ga)
        psi=1/(2*arb(m)).sqrt() if ga==1 else (arb(gg)/m).sqrt()
        for cnt,cc,mu in ((c['S'],fh/16,arb(1)),(c['I']+c['T'],fh,arb(1)/2)):
            al=2*PI*X0*(mu*m).sqrt()/ga
            MK2+=arb(1)/2*cnt*(2*arb(m)).sqrt()*cc*psi/arb(ga).sqrt()*(4*m*mu/D)**4*I(NU,al)
            F=PI**2*4*arb(m)**2*arb(ga)**2/3
            MK3+=arb(1)/2*cnt*cc*psi/(2*PI)/arb(ga)**2*(4*m*mu/D)**(arb(17)/4)*2*(mu/m).sqrt()*F*I(NU,al)
    # tails gamma>G for mock: I_nu(z/ga)<=ga^{-15/2} I_nu(z), cnt<=2ga, psi<=1
    al1=2*PI*X0*arb(m).sqrt()
    tailsum=arb(1)/(6*arb(G)**6)  # sum_{ga>G} ga^{-7} <= G^-6/6
    tail2=(2*arb(m)).sqrt()*fh*(1+arb(1)/16)*(4*arb(m)/D)**4*I(NU,al1)*tailsum*2
    tailsum3=arb(2)/(11*arb(G)**(arb(11)/2))  # sum_{ga>G} ga^{-13/2} <= G^{-11/2}/(11/2)
    tail3=fh*(1+arb(1)/16)/(2*PI)*(4*arb(m)/D)**(arb(17)/4)*2/arb(m).sqrt()*PI**2*4*arb(m)**2/3*I(NU,al1)*tailsum3*2
    MK2=(MK2+tail2)/ext; MK3=(MK3+tail3)/ext
    tot=E_main+E_miss+E_rs+E_TI+E2+MK2+MK3
    return dict(main=E_main,miss=E_miss,rs2=E_rs,TI=E_TI,E2=E2,MK2=MK2,MK3=MK3,total=tot)
PMIN=arb('0.538')  # pmin.py: min |P|^2 on |z|<=e^{-pi} >= 0.53816
if __name__=="__main__":
    import time
    args=sys.argv[1:]
    if args and args[0]=='-x':
        X0=arb(args[1]); args=args[2:]
    for m in [int(a) for a in args]:
        t=time.time(); r=certify(m)
        print(m,' '.join('%s=%s'%(k,v.str(3,radius=False)) for k,v in r.items()),'margin',(PMIN-r['total']).str(4,radius=False),'%.1fs'%(time.time()-t),flush=True)
