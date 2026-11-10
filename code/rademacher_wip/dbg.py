import sys
src=open("radmock2.py").read().rsplit("for n in (4,8,16,30):",1)[0]
exec(src)
# per-frame gamma=1,2 nonmock contributions
def radpart(n,l,GS):
    D=mp.mpf(2*m*n-l*l); out={}
    for ga in GS:
        for md in range(0,2*ga):
            de=-md
            if gcd(de,ga)!=1: continue
            al=pow(de% ga,-1,ga) if ga>1 else 0
            be=(al*de-1)//ga
            fr=frame((de,-be,-ga,al))
            data = pol['S'] if fr=='S' else pol['I']
            for (nt,lt,a) in data:
                if fr=='T': a=a*(-1)**int(2*nt)
                Dt=4*m*mp.mpf(nt.numerator)/nt.denominator - lt*lt
                ph=mp.expjpi(2*(al*Dt/(4*m)+de*(D/(4*m)))/ga)
                v=a*ph*weil(al,ga,de,lt,l)*(2*mp.pi/ga)*(-Dt/D)**(mp.mpf(15)/4)*mp.besseli(mp.mpf(15)/2,mp.pi*mp.sqrt(-Dt*D)/(m*ga))
                key=(ga,fr,str(nt),lt); out[key]=out.get(key,0)+v/2*(-1)**l*mp.sqrt(1j)
    return out
for n in (16,30):
  for l in (0,1):
    d=(-1)**(l+1)*coef(m,n,l)
    for k,v in sorted(radpart(n,l,(1,2,3)).items()): print(n,l,k,mp.nstr(v/d,6))
