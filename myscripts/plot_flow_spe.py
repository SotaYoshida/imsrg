#!/usr/local/bin/python3
import os, sys
import matplotlib.pyplot as plt
HOME=os.path.expanduser('~')
sys.path.append(HOME)
from python_myset.myplot import myplt

f = sys.argv[1]
fout = 'flow_spe_'+ os.path.splitext(os.path.basename(f))[0] + '.pdf'

plt.rcParams = myplt.set_style(roman=True)
fig,axs = myplt.set_canvas(1,1,6,8)
colors = myplt.colors()
proton_spe = {}
neutron_spe = {}

proton_s = {'0s1/2':3}
neutron_s = {'0s1/2':4}
proton_p = {'0p3/2':5, '0p1/2':7}
neutron_p = {'0p3/2':6, '0p1/2':8}
proton_sd = {'0d5/2':9, '0d3/2':11, '1s1/2':13}
neutron_sd = {'0d5/2':10, '0d3/2':12, '1s1/2':14}
proton_pf = {'0f7/2':15, '0f5/2':17, '1p3/2':19, '1p1/2':21}
neutron_pf = {'0f7/2':16, '0f5/2':18, '1p3/2':20, '1p1/2':22}
proton_sdg= {'0g9/2':23, '0g7/2':25, '1d5/2':27, '1d3/2':29, '2s1/2':31}
neutron_sdg= {'0g9/2':24, '0g7/2':26, '1d5/2':28, '1d3/2':30, '2s1/2':32}

proton_spe.update(proton_sd)
proton_spe.update(proton_pf)
neutron_spe.update(proton_sd)
neutron_spe.update(proton_pf)
cnt = 0
for key in neutron_spe.keys():
    label = key
    val = proton_spe[key]
    c = colors[val%len(colors)]
    X, Y = myplt.data2d(f,x=1,y=val)
    try:
        axs.plot(X[:-1],Y[:-1],label=label,c=c)
    except:
        pass
    cnt += 1


### comment ###
axs.set_xlabel('$s$')
axs.set_ylabel('Energy (MeV)')
#plt.annotate(r'$^{32}$Mg (g.s.; 0$^{+}$), NN+3N(400)', xy=(0,-25))
#plt.annotate(r'dash: neutron', xy=(38,-344))
#plt.annotate(r'solid: proton', xy=(38,-347))
### comment ###
plt.legend(ncol=3)
plt.tight_layout()
#plt.show()
plt.savefig(fout)
