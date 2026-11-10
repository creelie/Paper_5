# Figures for the CHL positivity paper (run from this directory; writes PNGs to ../)
import numpy as np, re, os, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
HERE=os.path.dirname(os.path.abspath(__file__)); OUT=os.path.join(HERE,'..')
CERT=os.path.join(HERE,'..','..','code','certificate')
plt.rcParams.update({'font.family':'serif','mathtext.fontset':'cm','font.size':10})
INK='#2b2b2b'; ACC='#c0392b'; BLUE='#2c6fbb'

def P(z,N=80):
    v=np.ones_like(z)
    for n in range(1,N+1): v=v*((1-z**n)*(1-z**(2*n)))**8
    return 1/v

# 1. |P(z)|^2 on the disc |z| <= e^{-pi}
r0=np.exp(-np.pi)
rr=np.linspace(0,r0,120); tt=np.linspace(0,2*np.pi,240)
R,T=np.meshgrid(rr,tt); Z=R*np.exp(1j*T); F=np.abs(P(Z))**2
fig=plt.figure(figsize=(6.4,4.6)); ax=fig.add_subplot(111,projection='3d')
ax.plot_surface(R*np.cos(T),R*np.sin(T),F,cmap=cm.viridis,linewidth=0,antialiased=True,alpha=0.92,rcount=120,ccount=240)
th=np.linspace(0,2*np.pi,400); zb=r0*np.exp(1j*th)
ax.plot(r0*np.cos(th),r0*np.sin(th),np.abs(P(zb))**2,color=INK,lw=1.2)
zmin=-r0; ax.scatter([-r0],[0],[abs(P(np.array([zmin+0j])))[0]**2],color=ACC,s=30,depthshade=False)
ax.text2D(0.04,0.93,r'minimum $|P|^2=0.5389$ at $z=-e^{-\pi}$ (red dot)',transform=ax.transAxes,color=ACC,fontsize=9)
ax.set_xlabel(r'$\mathrm{Re}\,z$',labelpad=6); ax.set_ylabel(r'$\mathrm{Im}\,z$',labelpad=6); ax.set_zlabel(r'$|P(z)|^2$',labelpad=4)
ax.set_zlim(0.4,2.3); ax.view_init(elev=24,azim=-58)
ax.set_xticks([-0.04,0,0.04]); ax.set_yticks([-0.04,0,0.04])
fig.tight_layout(); fig.savefig(os.path.join(OUT,'fig_pdisc.png'),dpi=220); plt.close(fig)

# 2. certificate: log10 of the error pieces, m = 3..100
rows=[]
for fn in ('cert_m3_7.out','cert_m8_100.out'):
    for line in open(os.path.join(CERT,fn)):
        m=int(line.split()[0]); d=dict(re.findall(r'(\w+)=([0-9.e+-]+)',line))
        rows.append((m,d))
rows.sort()
keys=[('MK3','mock terms'),('E2','$\\gamma\\geq2$'),('TI','frames $I,T$'),('rs2','walls $rs\\geq2$'),('miss','missing model terms'),('main','leading family')]
from mpl_toolkits.mplot3d.art3d import PolyCollection
fig=plt.figure(figsize=(6.8,4.8)); ax=fig.add_subplot(111,projection='3d')
cols=plt.cm.plasma(np.linspace(0.85,0.05,len(keys)))
verts=[]
for j,(k,lab) in enumerate(keys):
    ms=np.array([m for m,d in rows]); v=np.array([float(d[k])+(float(d['MK2']) if k=='MK3' else 0) for m,d in rows])
    lv=np.maximum(np.log10(np.maximum(v,1e-300)),-60)
    verts.append([(ms[0],-60)]+list(zip(ms,lv))+[(ms[-1],-60)])
pc=PolyCollection(verts,facecolors=cols,edgecolors=[INK]*len(keys),linewidths=0.6,alpha=0.8)
ax.add_collection3d(pc,zs=range(len(keys)),zdir='y')
ax.set_xlim(3,100); ax.set_ylim(-0.5,len(keys)-0.5)
ax.set_yticks(range(len(keys))); ax.set_yticklabels([l for k,l in keys],fontsize=7,rotation=-15,va='center',ha='left')
ax.set_xlabel(r'$Q^2=m$',labelpad=6); ax.set_zlabel(r'$\log_{10}$ (bound / extreme term)',labelpad=6)
ax.set_zlim(-60,0); ax.view_init(elev=22,azim=-62)
fig.tight_layout(); fig.savefig(os.path.join(OUT,'fig_certificate.png'),dpi=220); plt.close(fig)

# 3. margins
ms=np.array([m for m,d in rows]); tot=np.array([float(d['total']) for m,d in rows])
fig,ax=plt.subplots(figsize=(6.0,3.0))
ax.axhline(0.53816,color=INK,lw=1,ls='--'); ax.text(52,0.565,r'$\min_{|z|\leq e^{-\pi}}|P(z)|^2\geq0.53816$',fontsize=8,color=INK)
ax.fill_between(ms,0,tot,color=BLUE,alpha=0.35,lw=0); ax.plot(ms,tot,color=BLUE,lw=1.4,label='total error bound')
ax.plot(ms,0.53816-tot,color=ACC,lw=1.4,label='certified margin')
ax.set_xlim(3,100); ax.set_ylim(0,0.62); ax.set_xlabel(r'$Q^2=m$'); ax.legend(frameon=False,fontsize=8,loc='center right')
for s in ('top','right'): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(os.path.join(OUT,'fig_margins.png'),dpi=220); plt.close(fig)

# 4. schematic of walls: Farey geodesics q/s -> p/r in the omega upper half-plane
fig,ax=plt.subplots(figsize=(6.4,3.1))
from math import gcd
def arc(a,b,**kw):
    c=(a+b)/2; rad=abs(b-a)/2; t=np.linspace(0,np.pi,200); ax.plot(c+rad*np.cos(t),rad*np.sin(t),**kw)
for s in range(1,9):
    for r in range(1,9):
        if gcd(r,s)!=1: continue
        for q in range(-s,2*s+1):
            if (1+q*r)%s: continue
            p=(1+q*r)//s
            a,b=q/s,p/r
            if -0.05<=a and b<=1.05 and a<b:
                kw=dict(color=ACC,lw=1.8) if (r==1 and s==1) else dict(color='#7f8c8d',lw=0.6 if r*s>2 else 1.0)
                arc(a,b,**kw)
w0=0.37; dl=0.16
ax.plot([w0],[dl],'o',color=BLUE,ms=5); ax.text(w0+0.015,dl+0.012,r'$\omega_0+i\delta$',color=BLUE,fontsize=9)
ax.annotate('',xy=(w0,0.0),xytext=(w0,dl),arrowprops=dict(arrowstyle='-',color=BLUE,ls=':'))
ax.text(0.5,0.52,r'$r=s=1$: $G=(1\ 0;\ 1\ 1)$, leading family',color=ACC,ha='center',fontsize=8)
ax.text(0.72,0.2,r'$rs\geq2$: $\sqrt{t}\leq m/2+5$',color='#5f6a6a',fontsize=8)
ax.set_xlim(-0.03,1.03); ax.set_ylim(0,0.58); ax.set_aspect('equal'); ax.set_xlabel(r'$\omega=\tilde\ell/2m$')
ax.set_yticks([]); [ax.spines[s].set_visible(False) for s in ('top','right','left')]
fig.tight_layout(); fig.savefig(os.path.join(OUT,'fig_walls.png'),dpi=220); plt.close(fig)
print('ok')
