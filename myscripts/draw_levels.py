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
betalist=['0.5','1.0','1.1','1.2',"1.3","1.4","1.5"]
deltalist = [20]
x = 0
labels = []
for beta in betalist:
    for delta in deltalist:
        if(beta == '0.0'):
            snt = 'magnus_C12_custom_HF_hw25_e6_A12_delta10.snt'
        else:
            snt = 'magnus_C12_custom_HF_hw25_e6_A12_beta'+beta+'_delta10.snt'
        labels.append(r'$\beta=$'+beta)
        smf = 'summary_C12_'+os.path.splitext(os.path.basename(snt))[0]+'.txt'
        edict = plotlevel.get_energies_dct(smf,absolute=False,snt=snt)
        keys = []
        for key in edict.keys():
            if(key[2] > 1): keys.append(key)
        for key in keys:
            del edict[key]
        plotlevel.draw_energies(axs,edict,x,width,colorvar)
        if(x != 0):
            plotlevel.draw_connections(axs,edict_old,edict,x-1+width,x-width,colorvar)
        edict_old = edict
        x+=1
N = x - 1
plotlevel.put_JP_auto(axs,edict,N, 2, 0.15, colorvar)
plotlevel.set_frame(axs,range(x),labels)
plt.tight_layout()
plt.savefig('draw_levels.pdf')

