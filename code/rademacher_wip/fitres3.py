import sys
src=open("radmock3.py").read().rsplit("for n in (4,8,16,30):",1)[0]
exec(src)
prev=None
for n in (16,30,42):
    for l in (0,1):
        d=(-1)**(l+1)*coef(m,n,l)
        res=(d-rad(n,l)-mock(n,l,GM2)).real
        D=2*m*n-l*l
        print(n,l,D,mp.nstr(res,8), mp.nstr(mp.log(abs(res))/mp.sqrt(D),6))
