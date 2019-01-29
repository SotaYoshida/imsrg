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
e = 6
hw = 25
plt.rcParams = myplt.set_style()
fig,axs = myplt.set_canvas(1,1,10,9)

path = HOME+'path/to/file'
snt = 'snt file'
smf = 'kshell summary file'
sntf = path+snt
smff = path+smf
edict1 = {}
if(os.path.isfile(sntf) and os.path.isfile(smff)):
    edict1 = plotlevel.get_energies_dct(smff,absolute=True,snt=sntf)

path = HOME+'/path/to/file'
snt = 'snt file'
smf = 'kshell summary file'
sntf = path+snt
smff = path+smf
edict2 = {}
if(os.path.isfile(sntf) and os.path.isfile(smff)):
    edict2 = plotlevel.get_energies_dct(smff,absolute=True,snt=sntf)

plotlevel.draw_energies(axs,edict1,0,width,colorvar)
plotlevel.draw_energies(axs,edict2,1,width,colorvar)
plotlevel.draw_connections(axs,edict1,edict2,width,1-width,colorvar)

plotlevel.put_JP_auto(axs,edict2,1.4, 2, 0.15, colorvar)
plotlevel.set_frame(axs,[0,1],['$p$-shell','$psd$-shell'])

axs.set_xlim(-0.5,3)
plt.tight_layout()
plt.savefig('draw_levels.pdf')

