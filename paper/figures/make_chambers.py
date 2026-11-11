# Figure fig_chambers.png: the walls (even Farey edges) and chambers in the upper half-plane of tau_Y,
# the chamber K_0 and the vertical path of the proof of Theorem thm:Pinf (run from any directory; writes ../)
import os, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from math import gcd
HERE=os.path.dirname(os.path.abspath(__file__)); OUT=os.path.join(HERE,'..')
plt.rcParams.update({'font.family':'serif','mathtext.fontset':'cm','font.size':10})
INK='#2b2b2b'; ACC='#c0392b'; BLUE='#2c6fbb'; GREY='#7f8c8d'

X0,X1,YT,H=-0.12,1.12,1.15,0.32      # window in tau and height of the wall curtains
DEN=12
edges=[]                               # Farey edges a/c -- b/d (bd-ac... adjacent), finite ones
for c in range(1,DEN+1):
    for d in range(1,DEN+1):
        for a in range(-c,2*c+1):
            if gcd(a,c)!=1: continue
            for b in range(-d,2*d+1):
                if gcd(b,d)!=1 or b*c-a*d!=1: continue
                if a/c<X0-0.2 or b/d>X1+0.2: continue
                edges.append((a/c,b/d,(c*d)%2==0))
def arcpts(u,w,n=160):
    t=np.linspace(0,np.pi,n); ctr=(u+w)/2; r=(w-u)/2
    return ctr-r*np.cos(t), r*np.sin(t)

fig=plt.figure(figsize=(7.0,4.9)); ax=fig.add_subplot(111,projection='3d',computed_zorder=False)
# floor
ax.add_collection3d(Poly3DCollection([[(X0,0,0),(X1,0,0),(X1,YT,0),(X0,YT,0)]],facecolor='#f4f1ea',edgecolor='none',alpha=1.0,zorder=0))
# the chamber K_0: 0<x<1 outside the half-discs on (0,1/2) and (1/2,1)
xs=np.linspace(0,1,241)
low=np.where(xs<0.5, np.sqrt(np.maximum(0.0625-(xs-0.25)**2,0)), np.sqrt(np.maximum(0.0625-(xs-0.75)**2,0)))
poly=[(x,y,0.001) for x,y in zip(xs,low)]+[(1,YT,0.001),(0,YT,0.001)]
ax.add_collection3d(Poly3DCollection([poly],facecolor='#f5cba7',edgecolor='none',alpha=0.95,zorder=1))
ax.text(0.62,0.80,0.004,r'$K_0$',color=ACC,fontsize=15,ha='center',zorder=7)
# odd edges: dashed on the floor; even edges: curtains
for u,w,even in edges:
    if u<X0 or w>X1: continue
    x,y=arcpts(u,w)
    m=(y<=YT)
    if not m.any(): continue
    if not even:
        ax.plot(x[m],y[m],0*x[m],color=GREY,lw=0.6,ls=(0,(2,2)),zorder=2)
        continue
    hh=H*min(1.0,0.35+2.2*(w-u))
    verts=[(x[i],y[i],0) for i in np.where(m)[0]]+[(x[i],y[i],hh) for i in np.where(m)[0][::-1]]
    ax.add_collection3d(Poly3DCollection([verts],facecolor=BLUE,edgecolor='none',alpha=0.16 if w-u<0.2 else 0.28,zorder=5))
    ax.plot(x[m],y[m],hh+0*x[m],color=BLUE,lw=0.5 if w-u<0.2 else 1.0,zorder=6)
    ax.plot(x[m],y[m],0*x[m],color=BLUE,lw=0.7 if w-u<0.2 else 1.4,zorder=3)
# vertical walls Re tau = n (Farey edges n -- infinity)
for n in (0,1):
    ax.add_collection3d(Poly3DCollection([[(n,0,0),(n,YT,0),(n,YT,H),(n,0,H)]],facecolor=BLUE,edgecolor='none',alpha=0.22,zorder=5))
    ax.plot([n,n],[0,YT],[H,H],color=BLUE,lw=1.0,zorder=6); ax.plot([n,n],[0,YT],[0,0],color=BLUE,lw=1.4,zorder=3)
# the vertical path kappa* + i y and its crossings with walls
k=5/16
ax.plot([k,k],[YT,0.0],[0.006,0.006],color=ACC,lw=2.0,zorder=4)
ax.quiver(k,0.62,0.006,0,-0.14,0,color=ACC,arrow_length_ratio=0.5,lw=1.8,zorder=4)
hits=[]
for u,w,even in edges:
    if even and u<k<w:
        r=(w-u)/2; c=(u+w)/2; hits.append(np.sqrt(r*r-(k-c)**2))
for y in hits:
    ax.scatter([k],[y],[0.008],color=ACC,s=14,depthshade=False,zorder=4)
ax.text(k+0.03,1.0,0.01,r'$\kappa^*+i\mathrm{y}$',color=ACC,fontsize=10,zorder=7)
for v,lab in ((0,'0'),(0.5,r'$1/2$'),(1,'1')):
    ax.text(v,-0.07,0,lab,ha='center',va='top',fontsize=10,color=INK)
ax.set_xlim(X0,X1); ax.set_ylim(0,YT); ax.set_zlim(0,0.6)
ax.set_box_aspect((1.25,1.15,0.55))
ax.view_init(elev=34,azim=-64)
ax.set_xlabel(r'$\mathrm{Re}\,\tau_Y$',labelpad=4); ax.set_ylabel(r'$\mathrm{Im}\,\tau_Y$',labelpad=4)
ax.set_zticks([]); ax.set_xticks([0,0.5,1]); ax.set_xticklabels(['','','']); ax.set_yticks([0,0.5,1])
ax.xaxis.pane.set_visible(False); ax.yaxis.pane.set_visible(False); ax.zaxis.pane.set_visible(False)
ax.grid(False); ax.zaxis.line.set_visible(False)
fig.tight_layout(); fig.savefig(os.path.join(OUT,'fig_chambers.png'),dpi=220,bbox_inches='tight'); plt.close(fig)
from PIL import Image, ImageChops               # trim the empty margin left by the 3D axes
im=Image.open(os.path.join(OUT,'fig_chambers.png')).convert('RGB')
bb=ImageChops.difference(im,Image.new('RGB',im.size,(255,255,255))).getbbox()
im.crop((max(bb[0]-12,0),max(bb[1]-12,0),min(bb[2]+12,im.size[0]),min(bb[3]+12,im.size[1]))).save(os.path.join(OUT,'fig_chambers.png'),dpi=(220,220))
print('ok', len(edges), 'edges,', len(hits), 'walls crossed by the path')
