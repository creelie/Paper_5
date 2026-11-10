import sys
src=open("radmock4.py").read().rsplit("for n in (4,8,16,30):",1)[0]
exec(src)
for n in (10,20):
    for l in range(m+1):
        d=(-1)**(l+1)*coef(m,n,l); r=rad(n,l); mk=mock(n,l,GM2)
        print(n,l,'d',d,'rad/d',mp.nstr(r.real/d,10),'mock/d',mp.nstr(mk.real/d,6),'res/d',mp.nstr((d-r-mk).real/d,6))
