#!/usr/local/bin/python3
import os, sys
import numpy as np
import matplotlib.pyplot as plt
HOME=os.path.expanduser('~')
sys.path.append(HOME)
from python_myset.myplot import myplt
import matplotlib.animation as animation
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

fe = sys.argv[1]
fm = sys.argv[2]
fout = os.path.splitext(os.path.basename(fe))[0] + '.mp4'
plt.rcParams = myplt.set_style(roman=True)
fig,axs = myplt.set_canvas_grid([(10,4),(5,4),(5,4)],width=10,height=8,wspace=5,hspace=4)
X,Y=myplt.data2d(fe,x=1,y=2)
data = np.loadtxt(fm)

axs[0].plot(X,Y)
dshape = data.shape
Nmat = dshape[0] // dshape[1]
Mats = data.reshape(Nmat//2,2,dshape[1],dshape[1])
ims=[]
for i in range(Nmat//2-1):
    imp, = axs[0].plot(X[i],Y[i],ms=10,c='r',marker='o')
    imv = axs[1].imshow(Mats[i,0,:,:], animated=True, vmin=-1, vmax=1, cmap='seismic')
    ime = axs[2].imshow(Mats[i,1,:,:], animated=True,
            norm=colors.SymLogNorm(linthresh=1.e-6,linscale=1,vmin=-1.e-3,vmax=1.e-3), cmap='seismic')
    ims.append([imp,imv,ime])
ani = animation.ArtistAnimation(fig, ims, interval=100, blit=True,repeat_delay=1000)
axs[0].set_xlabel('s')
axs[0].set_ylabel('Energy (MeV)')
axs[1].set_xlabel('$\Gamma$')
axs[2].set_xlabel('$\eta$')
divider1=make_axes_locatable(axs[1])
divider2=make_axes_locatable(axs[2])
cax1 = divider1.append_axes('right','5%',pad='3%')
cax2 = divider2.append_axes('right','5%',pad='3%')
fig.colorbar(imv,ax=axs[1], cax = cax1)
fig.colorbar(ime,ax=axs[2], cax = cax2)
axs[1].margins(2,2)
plt.tight_layout()
ani.save(fout)
print("See output: " + fout)





