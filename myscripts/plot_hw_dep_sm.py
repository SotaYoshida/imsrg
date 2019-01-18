#!/usr/local/bin/python3
import sys, os
import matplotlib.pyplot as plt
HOME=os.path.expanduser('~')
sys.path.append(HOME)
from python_myset.myplot import myplt
from python_myset.SMint.kshell import kshell

colors = myplt.colors()
marks = myplt.markers()
elist = [2,4,6,8]
hwlist = [10,15,20,25,30,35,40,45,50]
cnt = 0
for e in elist:
    c = colors[cnt]
    m = marks[cnt]
    x = []
    y = []
    for hw in hwlist:
        snt = 'magnus_O18_custom_NAT_hw' + str(hw) +\
                '_e' + str(e) + '_A18.snt'
        smf = 'summary_O18_magnus_O18_custom_NAT_hw' + str(hw) +\
                '_e' + str(e) + '_A18.txt'
        if(not os.path.isfile(snt)): continue
        if(not os.path.isfile(smf)): continue
        k = kshell()
        k.read(snt)
        lvl = myplt.data_point(smf,r=5,c=5) # ground-state
        x.append(hw)
        y.append(lvl + k.zero_body)
    plt.plot(x,y,c=c,marker=m)
    cnt += 1
plt.show()

