import sys
exec(open('radfull.py').read().split("for n in (4,8,16,30):")[0].replace("m=int(sys.argv[1]); GMAX=int(sys.argv[2])","m=int(sys.argv[1]); GMAX=int(sys.argv[2]); PW=float(sys.argv[3]) if len(sys.argv)>3 else 0"))
def fhat(m): return D1(m)
def kl(D,mu_t,ga,lt,l,frames):
    # sum over -delta in [0,2ga), frame of Gamma^{-1} in frames; returns dict frame->sum
    out={}
    for md in range(0,2*ga):
        de=-md
        if gcd(de,ga)!=1: continue
        al=pow(de%ga,-1,ga) if ga>1 else 0
        be=(al*de-1)//ga
        fr=frame((de,-be,-ga,al))
        mu,c=frames[fr]
        mt=mu_t(mu)
        ph=mp.expjpi(2*(al*mt+de*D/(4*m))/ga)
        out[fr]=out.get(fr,0)+c*ph*weil(al,ga,de,lt,l)
    return out
def f_(ga,g,x):
    if g%(2*m*ga)==0: return mp.pi**2/mp.sinh(mp.pi*x/ga)**2 - mp.mpf(ga)**2/x**2
    return mp.pi**2/mp.sinh(mp.pi*x/ga - 1j*mp.pi*g/(2*m*ga))**2
def mock(n,l,GM):
    D=mp.mpf(2*m*n-l*l); fm=fhat(m)
    frames={'S':(mp.mpf(-1),mp.mpf(fm)/16*SW),'I':(mp.mpf(-1)/2,mp.mpf(fm)*IW),'T':(mp.mpf(-1)/2,-mp.mpf(fm)*IW)}
    tot=0
    for ga in range(1,GM+1):
        # line 2
        for fr,(mu,c) in frames.items():
            pass
        K=kl(D,lambda mu:mu,ga,0,l,frames)
        for fr,v in K.items():
            mu=-frames[fr][0]
            tot+= -mp.sqrt(2*m)*v/mp.sqrt(ga)*(4*m*mu/D)**4*mp.besseli(8,4*mp.pi*mp.sqrt(mu*D/(4*m))/ga)*mu**PW
        # line 3
        for g in range(-m*ga,m*ga):
            lt=g%(2*m)
            for fr0 in ('S','I','T'):
                sub={k:(v if k==fr0 else (v[0],0)) for k,v in frames.items()}
                K=kl(D,lambda mu:mu-mp.mpf(g*g)/(4*m),ga,lt,l,sub)
                v=K.get(fr0,0)
                if v==0: continue
                mu=-frames[fr0][0]
                X=mp.sqrt(mu/m)
                integ=mp.quad(lambda x: f_(ga,g,x)*(1-m*x*x/mu)**(mp.mpf(17)/4)*mp.besseli(mp.mpf(17)/2,4*mp.pi*mp.sqrt((mu-m*x*x)*D/(4*m))/ga),[-X,0,X])
                tot+= 1/(2*mp.pi)*v/ga**2*(4*m*mu/D)**(mp.mpf(17)/4)*integ*mu**PW
    return tot/2*(-1)**l*mp.sqrt(1j)
GM2=int(sys.argv[4]); SW=float(sys.argv[5]); IW=float(sys.argv[6])
for n in (4,8,16,30):
    for l in range(0,m+1):
        d=(-1)**(l+1)*coef(m,n,l)
        r=rad(n,l); mk=mock(n,l,GM2)
        print(n,l,'diff/d',mp.nstr((d-r)/d,6),'mock/d',mp.nstr(mk/d,6),'(d-r-mk)/d',mp.nstr((d-r-mk)/d,6))
