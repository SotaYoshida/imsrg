#!/usr/local/bin/python3
import sys, os
import matplotlib.pyplot as plt
HOME=os.path.expanduser('~')
sys.path.append(HOME)
from python_myset.myplot import myplt
from python_myset.myplot import plotlevel
from python_myset.myplot.SMint.kshell import kshell

width=0.2
colorvar = None
plt.rcParams = myplt.set_style()
fig,axs = myplt.set_canvas(1,1,width=20,height=8)
betalist=['0.0','0.2','0.4','1.0']
deltalist = range(15,40,5)
x = 0
labels = []
for beta in betalist:
    for delta in deltalist:
        path = 'betacm0/'
        if(beta == '0.0'):
            snt = 'magnus_O16_custom_HF_hw25_e8_A16_delta'+str(delta)+'.snt'
        else:
            snt = 'magnus_O16_custom_HF_hw25_e8_A16_beta'+str(beta)+\
                    '_delta'+str(delta)+'.snt'
        labels.append(r'$\beta=$'+beta+r', $\Delta=$'+str(delta))
        smf = path+'summary_O16_'+os.path.splitext(os.path.basename(snt))[0]+'.txt'
        edict = plotlevel.get_energies_dct(smf,absolute=True,snt=snt)
        plotlevel.draw_energies(axs,edict,x,width,colorvar)
        if(x != 0):
            plotlevel.draw_connections(axs,edict_old,edict,x-1+width,x-width,colorvar)
        edict_old = edict
        x+=1
N = x - 1
plotlevel.put_JP_auto(axs,edict,N, 2, 0.15, colorvar)
plotlevel.set_frame(axs,range(x),labels)
plt.tight_layout()
plt.savefig('draw_levels_betacm0.pdf')

